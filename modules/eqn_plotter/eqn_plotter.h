#ifndef __EQN_PLOTTER_H__
#define __EQN_PLOTTER_H__

#include <cmath>
#include <fstream>

#include "mio/altro/types.h"
#include "mio/math/math.h"
#include "mio/altro/algorithm.h"
#include "mio/math/ransac.h"
#include "mio/math/stats.h"
#include "mio/qt/qwt_graph.h"

#include "cook_torrance.h"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"


class CPlotBeckmann : public CMySliderGraph{
  public:
    CPlotBeckmann(){
      SetNumSlider(1); //mu
      SetupSlider(0, 0, 1.0, 0.5, 3, "mu"); //mu
      graph->qwt_plot->setTitle("Beckmann Distribution");
      graph->SetNumCurve(1);
      graph->SetCurveTitle(0, "B()");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_step = 500;
      std::vector<double> phi_vec(num_step), y_data_vec(num_step);
      mio::FillRange(phi_vec, -M_PI/2, M_PI/2);
      const double mu = SliderValue(0);

      for(size_t i = 0; i < num_step; ++i)
        y_data_vec[i] = BeckmannDistribution(mu, phi_vec[i]);

      graph->SetCurveData(0, phi_vec.data(), y_data_vec.data(), num_step);
    }
};


class CPlotFresnelDerivatives : public CMySliderGraph{
  public:
    CPlotFresnelDerivatives(){
      SetNumSlider(1); //theta
      SetupSlider(0, 0, M_PI/2.0, M_PI/4.0, 3, "theta"); //theta
      graph->qwt_plot->setTitle("Fresnel n derivative");
      graph->SetNumCurve(8);
      graph->SetCurveTitle(0, "F(n)_sim");
      graph->SetCurveTitle(1, "F(n)_app");
      graph->SetCurveTitle(2, "F'(n)_sim");
      graph->SetCurveTitle(3, "F'(n)_app");
      graph->SetCurveTitle(4, "F(k)_sim");
      graph->SetCurveTitle(5, "F(k)_app");
      graph->SetCurveTitle(6, "F'(k)_sim");
      graph->SetCurveTitle(7, "F'(k)_app");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_step = 500;
      //n_vec and k_vec are x-axis data. the others are y-axis data.
      std::vector<double> n_vec(num_step), k_vec(num_step),
                          FoN_app_vec(num_step), FoN_app_prime_vec(num_step),
                          FoN_sim_vec(num_step), FoN_sim_prime_vec(num_step),
                          FoK_app_vec(num_step), FoK_app_prime_vec(num_step),
                          FoK_sim_vec(num_step), FoK_sim_prime_vec(num_step);
      std::vector<double> FoN_complex_vec(num_step), FoK_complex_vec(num_step);

      K_DBL min_n = 0.25, max_n = 10.0;
      K_DBL n_step = (max_n - min_n) / static_cast<double>(num_step);
      K_DBL min_k = 0.25, max_k = 10.0;
      K_DBL k_step = (max_k - min_k) / static_cast<double>(num_step);
      K_DBL err_a = 1.0, err_gamma = 1.0; //cos(theta = 0.0) = 1
      double n = min_n,
             k = min_k,
             d_dn, d_dk, temp;
      const double theta = SliderValue(0);

      for(size_t i = 0; i < num_step; ++i){
        n_vec[i] = n;
        k_vec[i] = k;
        // Fresnel Simple
        FoN_sim_vec[i] = FresnelSimple(n, 5.0, theta);
        FoK_sim_vec[i] = FresnelSimple(1.5, k, theta);
        FresnelSimple_deriv(n, 5.0, theta, d_dn, temp);
        FresnelSimple_deriv(1.5, k, theta, temp, d_dk);
        FoN_sim_prime_vec[i] = d_dn;
        FoK_sim_prime_vec[i] = d_dk;
        // Fresnel Approximation
        FoN_app_vec[i] = FresnelMetalApprox(n, 5.0, theta, err_a, err_gamma);
        FoK_app_vec[i] = FresnelMetalApprox(1.5, k, theta, err_a, err_gamma);
        FresnelMetalApprox_deriv(n, 5.0, theta, d_dn, temp);
        FresnelMetalApprox_deriv(1.5, k, theta, temp, d_dk);
        FoN_app_prime_vec[i] = d_dn;
        FoK_app_prime_vec[i] = d_dk;
        n += n_step;
        k += k_step;
      }

      const int flags = MyGraph::DEFAULT_RAND_COLOR | MyGraph::COPY_CURVE_DATA;// | MyGraph::FIT_SPLINE;
      graph->SetCurveData(0, n_vec.data(), FoN_sim_vec.data(), num_step, flags);
      graph->SetCurveData(1, n_vec.data(), FoN_app_vec.data(), num_step, flags);
      graph->SetCurveData(2, n_vec.data(), FoN_sim_prime_vec.data(), num_step, flags);
      graph->SetCurveData(3, n_vec.data(), FoN_app_prime_vec.data(), num_step, flags);
      graph->SetCurveData(4, k_vec.data(), FoK_sim_vec.data(), num_step, flags);
      graph->SetCurveData(5, k_vec.data(), FoK_app_vec.data(), num_step, flags);
      graph->SetCurveData(6, k_vec.data(), FoK_sim_prime_vec.data(), num_step, flags);
      graph->SetCurveData(7, k_vec.data(), FoK_app_prime_vec.data(), num_step, flags);
    }
};


class CPlotFresnel : public CMySliderGraph{
  public:
    CPlotFresnel(){
      SetNumSlider(2); //n, k
      SetupSlider(0, 0.25, 10.0, 1.5, 3, "n"); //n
      SetupSlider(1, 0.25, 10.0, 5.0, 3, "k"); //k
      graph->qwt_plot->setTitle("Fresnel");
      graph->SetNumCurve(3);
      graph->SetCurveTitle(0, "F()_cmplx");
      graph->SetCurveTitle(1, "F()_simpl");
      graph->SetCurveTitle(2, "F()_apprx");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_step = 500;
      std::vector<double> theta_vec(num_step), fresnel_complex(num_step),
                          fresnel_simple(num_step), fresnel_approx(num_step);
      mio::FillRange(theta_vec, 0.0, M_PI/2.0);
      const double n = SliderValue(0),
                   k = SliderValue(1);

      for(size_t i = 0; i < num_step; ++i){
        fresnel_complex[i] = FresnelComplex(n, k, theta_vec[i]);
        fresnel_simple[i] = FresnelSimple(n, k, theta_vec[i]);
        fresnel_approx[i] = FresnelMetalApprox(n, k, theta_vec[i], 2.0*n, 7.0);
      }

      graph->SetCurveData(0, theta_vec.data(), fresnel_complex.data(), num_step);
      graph->SetCurveData(1, theta_vec.data(), fresnel_simple.data(), num_step);
      graph->SetCurveData(2, theta_vec.data(), fresnel_approx.data(), num_step);
    }
};


class CPlotCookTorrance : public CMySliderGraph{
  struct CookTorranceData{
    vec3d_t N, V, H, L;
  };

  public:
    CPlotCookTorrance(){
      SetNumSlider(3); //n, k, mu
      SetupSlider(0, 0.25, 10.0, 1.5, 3, "n"); //n
      SetupSlider(1, 0.25, 10.0, 5.0, 3, "k"); //k
      SetupSlider(2, 0, 1.0, 0.5, 3, "mu"); //mu
      graph->qwt_plot->setTitle("Cook Torrance");
      graph->SetNumCurve(1);
      graph->SetCurveTitle(0, "R()");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_train_sample = 500;
      std::vector<CookTorranceData> ctd_vec(num_train_sample);
      std::vector<double> y_data_vec(num_train_sample), theta_vec(num_train_sample);
      const vec3d_t V(0.0, 0.0, 1.0),
                    L( 0.0, std::sin( sm::DegToRad(10) ), std::cos( sm::DegToRad(10) ) );
      const vec3d_t H = sm::HalfwayVectorNorm3(L, V);
      mio::FillRange(theta_vec, sm::DegToRad(0.0), sm::DegToRad(90.0));
      const double n = SliderValue(0),
                   k = SliderValue(1),
                   mu = SliderValue(2);

      for(size_t i = 0; i < num_train_sample; ++i){
        //simulate flat plate that is rotation, camera and illuminator are stationary
        ctd_vec[i].N = vec3d_t( 0.0, std::sin(theta_vec[i]), std::cos(theta_vec[i]) );
        ctd_vec[i].V = V;
        ctd_vec[i].L = L;
        ctd_vec[i].H = H;
      }
      for(size_t i = 0; i < num_train_sample; ++i)
        y_data_vec[i] = CookTorrance(ctd_vec[i].N, ctd_vec[i].V, ctd_vec[i].H, ctd_vec[i].L, n, k, mu);

      graph->SetCurveData(0, theta_vec.data(), y_data_vec.data(), num_train_sample);
    }
};


class CPlotPlanck : public CMySliderGraph{
  const double  h = 6.6260695729e-34, //Planck's constant 6.62606957(13) × 10-34 m^2 kg / s
               kB = 1.380648813e-23,  //Boltzmann constant 1.3806488(29) × 10-23 m^2 kg s^-2 K^-1
                c = 299792458;

  double SolvePlanckLaw(double lambda, double temperature){
    return (2.0f*h*c*c / pow(lambda, 5)) / (exp( (h*c) / (lambda*kB*temperature) ) - 1.0f);
  }

  public:
    CPlotPlanck(){
      SetNumSlider(1);
      SetupSlider(0, 525, 775, 650, 3, "K"); //T
      graph->qwt_plot->setTitle("Planck Black Body (850nm to 1700nm)");
      graph->SetNumCurve(1);
      graph->SetCurveTitle(0, "Planck");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_train_sample = 500;
      std::vector<double> y_data_vec(num_train_sample), lambda_vec_m(num_train_sample), lambda_vec_nm(num_train_sample);
      mio::FillRange<double>(lambda_vec_m, 850e-9, 1700e-9);
      mio::FillRange<double>(lambda_vec_nm, 850, 1700);
      const double T = SliderValue(0);

      for(size_t i = 0; i < num_train_sample; ++i)
        y_data_vec[i] = SolvePlanckLaw(lambda_vec_m[i], T);

      graph->SetCurveData(0, lambda_vec_nm.data(), y_data_vec.data(), num_train_sample);
    }
};


class CPlotSuperPos : public CMySliderGraph{
  public:
    CPlotSuperPos(){
      SetNumSlider(4);
      SetupSlider(0, 0.0, 180.0, 0.0, 1, "phase_2"); //phase
      SetupSlider(1, 0.0, 2.0, 1.0, 2, "amp_2"); //amplitude
      SetupSlider(2, 0.25, 10.0, 1.0, 2, "freq_1"); //frequency wave 1
      SetupSlider(3, 0.25, 10.0, 1.0, 2, "freq_2"); //frequency wave 2
      graph->qwt_plot->setTitle("Superposition");
      graph->SetNumCurve(3);
      graph->SetCurveTitle(0, "wave_1");
      graph->SetCurveTitle(1, "wave_2");
      graph->SetCurveTitle(2, "wave_1+2");

      InitUI();
    }

    void PlotCurves(){
      const size_t num_train_sample = 500;
      std::vector<double> theta_vec(num_train_sample), theta_deg_vec(num_train_sample),
                          wave_1(num_train_sample), wave_2(num_train_sample), wave_12(num_train_sample);
      mio::FillRange<double>(theta_vec, -6.0*M_PI, 6.0*M_PI);
      mio::FillRange<double>(theta_deg_vec, sm::RadToDeg(-4.0*M_PI), sm::RadToDeg(4.0*M_PI));
      const double phase = SliderValue(0),
                   amplitude = SliderValue(1),
                   frequency_1 = SliderValue(2),
                   frequency_2 = SliderValue(3);

      for(size_t i = 0; i < num_train_sample; ++i){
        wave_1[i] = std::sin(theta_vec[i]*frequency_1);
        wave_2[i] = std::sin((theta_vec[i]+sm::DegToRad(phase))*frequency_2)*amplitude;
        wave_12[i] = wave_1[i] + wave_2[i];
      }

      graph->SetCurveData(0, theta_deg_vec.data(), wave_1.data(), num_train_sample);
      graph->SetCurveData(1, theta_deg_vec.data(), wave_2.data(), num_train_sample);
      graph->SetCurveData(2, theta_deg_vec.data(), wave_12.data(), num_train_sample);
    }
};


class CFilterData : public CMySliderGraph{
  std::vector<double> y_data_raw_vec_, x_data_raw_vec_;

  public:
    CFilterData(const std::string kFileName){
      SetNumSlider(2);
      SetupSlider(0, 0.0, 1, 0.0050, 4, "max_gradient");
      SetupSlider(1, 0.0, 5, 0.850, 3, "std_dev_coef");
      graph->qwt_plot->setTitle("Merc Data Filt Test");
      graph->SetNumCurve(3);
      graph->SetCurveTitle(0, "raw");
      graph->SetCurveTitle(1, "grad_filter");
      graph->SetCurveTitle(2, "stat_filter");

      std::fstream in_file(kFileName, std::ios_base::in);
      float y;
      while(in_file >> y)
        y_data_raw_vec_.push_back(y);
      std::cout << "read " << y_data_raw_vec_.size() << " points\n";
      x_data_raw_vec_.resize(y_data_raw_vec_.size());
      std::iota(x_data_raw_vec_.begin(), x_data_raw_vec_.end(), 0);

      InitUI();
    }

    void GradientFilter(const std::vector<double> &src_data, std::vector<double> &filt_data,
                        std::vector<size_t> &filt_data_idx, const double max_gradient){
      filt_data.clear();
      filt_data.reserve(src_data.size()-2);
      const size_t kDataInSizeMinusOne = src_data.size()-1;
      for(size_t i = 0; i < kDataInSizeMinusOne; ++i){
        double gradient = 1 + src_data[i]*src_data[i+1];
        std::cout << gradient << " ";
      }
      std::cout << std::endl;
    }

    void PlotCurves(){
      std::vector<double> y_data_grad_filt_vec, x_data_grad_filt_vec,
                          y_data_stat_filt_vec, x_data_stat_filt_vec;
      std::vector<size_t> filt_idx;
      const double max_gradient = SliderValue(0),
                   std_dev_coef = SliderValue(1);

      sm::StatisticalOutlierRemoval(y_data_raw_vec_, y_data_stat_filt_vec, filt_idx, std_dev_coef);
      std::copy(filt_idx.begin(), filt_idx.end(), std::back_inserter(x_data_stat_filt_vec));
      GradientFilter(y_data_stat_filt_vec, y_data_grad_filt_vec, filt_idx, max_gradient);
      std::copy(filt_idx.begin(), filt_idx.end(), std::back_inserter(x_data_grad_filt_vec));

      //double sum = 0;
      //for(const auto &val : y_data_stat_filt_vec)
        //sum += val;
      //std::cout << "mean=" << sum/static_cast<double>(y_data_stat_filt_vec.size()) << std::endl;

      graph->SetCurveData(0, x_data_raw_vec_.data(), y_data_raw_vec_.data(), x_data_raw_vec_.size());
      graph->SetCurveData(1, x_data_grad_filt_vec.data(), y_data_grad_filt_vec.data(), x_data_grad_filt_vec.size());
      graph->SetCurveData(2, x_data_stat_filt_vec.data(), y_data_stat_filt_vec.data(), x_data_stat_filt_vec.size());
    }
};


class CPlotHesData : public CMySliderGraph{
  std::vector<double> y_vec_, x_vec_;
  std::vector<vertex2d_t> pnt_vec_;

  public:
    CPlotHesData(const std::string kFileName){
      SetNumSlider(1);
      // (slider_idx, min, max, value, num_dec = 3, label_text = QString())
      SetupSlider(0, 3, 151, 15, 0, "kernel size");
      graph->qwt_plot->setTitle("HES Data Filt Test");
      graph->SetNumCurve(4);
      graph->SetCurveTitle(0, "raw");
      graph->SetCurveTitle(1, "filt");
      graph->SetCurveTitle(2, "line_a");
      graph->SetCurveTitle(3, "line_b");

      std::fstream in_file(kFileName, std::ios_base::in);
      float x = 0, y;
      while(in_file >> y){
        y_vec_.push_back(y);
        x_vec_.push_back(x);
        pnt_vec_.push_back(vertex2d_t(x++, y));
      }
      std::cout << "read " << pnt_vec_.size() << " points\n";

      InitUI();
    }

    void PlotCurves(){
      std::vector<double> blur_vec;
      int kernel_width = static_cast<int>(SliderValue(0)) % 2 == 0 ? SliderValue(0)+1 : SliderValue(0);
      cv::blur(y_vec_, blur_vec, cv::Size(kernel_width, 1));

      std::vector< std::pair<uint32_t, line3d_t> > detected_lines;
      const double kMinInlierDist = 2;
      const uint32_t kMinNumInlier = 200;
      std::vector< std::vector<uint32_t> > inlier_index_vvec;
      
      std::vector<CVertex2d> blur_pnt_vec;
      CVertex2d::Merge(x_vec_, blur_vec, blur_pnt_vec);
      sm::RansacDetectLines2D(blur_pnt_vec, detected_lines, kMinInlierDist, kMinNumInlier, inlier_index_vvec);

      graph->SetCurveData(0, x_vec_.data(), y_vec_.data(), x_vec_.size(), MyGraph::DEFAULT_RAND_COLOR);
      graph->SetCurveData(1, x_vec_.data(), blur_vec.data(), x_vec_.size());

      std::vector<double> y_vec, x_vec;
      for(size_t i = 0; i < detected_lines.size() && (i+2) < 4; ++i){
        int min = INT_MAX, max = INT_MIN;
        for(size_t j = 0; j < inlier_index_vvec[i].size(); ++j){
          if(inlier_index_vvec[i][j] < min)
            min = inlier_index_vvec[i][j];
          if(static_cast<int>(inlier_index_vvec[i][j]) > max)
            max = inlier_index_vvec[i][j];
        }
        min -= 25;
        max += 25;
        x_vec.push_back(min);
        y_vec.push_back(sm::SolveLineEqn2(detected_lines[i].second, min));
        x_vec.push_back(max);
        y_vec.push_back(sm::SolveLineEqn2(detected_lines[i].second, max));
        graph->SetCurveData(i+2, x_vec.data(), y_vec.data(), x_vec.size());
        x_vec.clear();
        y_vec.clear();
      }
    }
};


class CPlotHesDataB : public CMySliderGraph{
  std::vector<double> y_vec_, x_vec_;
  std::vector<vertex2d_t> pnt_vec_;

  public:
    CPlotHesDataB(const std::string kFileName){
      SetNumSlider(5);
      // (slider_idx, min, max, value, num_dec = 3, label_text = QString())
      SetupSlider(0, 0, 1, 0.25, 2, "min_inlier_percent");
      SetupSlider(1, 0, 3000, 750, 0, "num_lead_sample");
      SetupSlider(2, 0, 1000, 200, 0, "num_tail_sample");
      SetupSlider(3, 0, 10, 2, 1, "std_dev_coef");
      SetupSlider(4, 0, 100, 50, 0, "update");
      graph->qwt_plot->setTitle("HES Data Filt Test");
      graph->SetNumCurve(5);
      graph->SetCurveTitle(0, "raw");
      graph->SetCurveTitle(1, "line_a");
      graph->SetCurveTitle(2, "line_b");
      graph->SetCurveTitle(3, "line_c");
      graph->SetCurveTitle(4, "line_d");

      std::fstream in_file(kFileName, std::ios_base::in);
      float x = 0, y;
      while(in_file >> y){
        y_vec_.push_back(y);
        x_vec_.push_back(x);
        pnt_vec_.push_back(vertex2d_t(x++, y));
      }
      std::cout << "read " << pnt_vec_.size() << " points\n";

      InitUI();
    }

    void PlotCurves(){
      graph->SetCurveData(0, x_vec_.data(), y_vec_.data(), x_vec_.size(), MyGraph::DEFAULT_RAND_COLOR);

      std::vector<double> lead_y_vec(y_vec_.begin(), y_vec_.begin()+static_cast<size_t>(SliderValue(1)));
      std::vector<vertex2d_t> lead_pnt_vec(pnt_vec_.begin(), pnt_vec_.begin()+static_cast<size_t>(SliderValue(1)));
      std::vector<double> tail_y_vec(y_vec_.end()-SliderValue(2), y_vec_.end());
      std::vector<vertex2d_t> tail_pnt_vec(pnt_vec_.end()-SliderValue(2), pnt_vec_.end());

      double min, max, variance, std_dev, mean, median;
      std::set<double> set_data;
      sm::ComputeStatistics(lead_y_vec, min, max, variance, std_dev, mean, median, set_data);
      //std:: cout << "variance: " << variance << ", std_dev: " << std_dev << std::endl;

      std::vector< std::pair<uint32_t, line3d_t> > detected_lines;
      const double kMinInlierDist = SliderValue(3)*std_dev;
      const uint32_t kMinNumInlier = SliderValue(0)*SliderValue(1);
      std::vector< std::vector<uint32_t> > inlier_index_vvec;
      sm::RansacDetectLines2D(lead_pnt_vec, detected_lines, kMinInlierDist, kMinNumInlier, inlier_index_vvec);

      // Plot up to four lines
      std::vector<double> y_vec, x_vec;
      for(size_t i = 0; i < detected_lines.size() && i < 4; ++i){
        x_vec.push_back(x_vec_.front());
        y_vec.push_back(sm::SolveLineEqn2(detected_lines[i].second, x_vec_.front()));
        x_vec.push_back(x_vec_.back());
        y_vec.push_back(sm::SolveLineEqn2(detected_lines[i].second, x_vec_.back()));
        graph->SetCurveData(i+1, x_vec.data(), y_vec.data(), x_vec.size());
        x_vec.clear();
        y_vec.clear();
      }
      for(size_t i = detected_lines.size()-1; i < 4; ++i)
        graph->SetCurveData(i+1, x_vec.data(), y_vec.data(), x_vec.size());
    }
};


/*
void SimpleOneBlock(const double weight, const double shaft_length,
                    const double el_coef, const double s_coef){
  std::cout << "Simple Supported, One Block: ";
  std::cout << ( ( weight * std::pow(shaft_length, 3) ) / (48.0 * el_coef) ) +
               ( ( 5.0 * s_coef * std::pow(shaft_length, 4) ) / (384.0 * el_coef) ) << std::endl;
}

void SimpleTwoBlock(const double weight, const double shaft_length, const double block_separation,
                    const double el_coef, const double s_coef){
  const double a = shaft_length*0.5 - block_separation;
  std::cout << "Simple Supported, Two Block: ";
  std::cout << ( ( weight * a * ( (3.0 * shaft_length * shaft_length) - (4.0 * a * a) ) ) / (48.0 * el_coef) ) +
               ( ( 5.0 * s_coef * std::pow(shaft_length, 4) ) / (384.0 * el_coef) ) << std::endl;
}
*/

#endif //__EQN_PLOTTER_H__

