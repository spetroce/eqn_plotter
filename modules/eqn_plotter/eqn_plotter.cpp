#include <qapplication.h>
#include <QStyleFactory>
#include "eqn_plotter.h"


int main(int argc, char *argv[]){
  QApplication a(argc, argv);
  QApplication::setStyle(QStyleFactory::create("Fusion"));
  STD_INVALID_ARG_E(argc >= 2)

  MyGraph *mg;
  CMySliderGraph *sg;

  if(std::atoi(argv[1]) < 0){
    std::vector<double> x(32), y(32);
    std::iota(x.begin(), x.end(), 0);
    const double step = M_PI/31.;
    double val = -step;
    std::generate(y.begin(), y.end(), [&]{ val+=step; return 250.*std::sin(val);});
    mg = new MyGraph();
    mg->SetNumCurve(1);
    mg->SetCurveData(0, x.data(), y.data(), 32);
    mg->SetPlotSize(640, 480);
    mg->qwt_plot->show();
  }
  else{
    switch(std::atoi(argv[1])){
      case 0:
        sg = new CPlotFresnel();
        break;
      case 1:
        sg = new CPlotFresnelDerivatives();
        break;
      case 2:
        sg = new CPlotBeckmann();
        break;
      case 3:
        sg = new CPlotCookTorrance();
        break;
      case 4:
        sg = new CPlotPlanck();
        break;
      case 5:
        sg = new CPlotSuperPos();
        break;
      case 6:
        sg = new CFilterData(std::string(argv[2]));
        break;
      case 7:
        sg = new CPlotHesData(std::string(argv[2]));
        break;
      default:
        printf("invalid plot type %d\n", std::atoi(argv[1]));
        return 0;
    }
    sg->show();
  }

  return a.exec();
}

