#include <iostream>
#include <functional>
#include <vector>
#include <random>
#include <mutex>
#include <chrono>
#include <thread>

#define GET_VARIABLE_NAME(Variable) (#Variable)

class ThreadJoiner
{
public:
    ThreadJoiner(std::thread t)
            :t(std::move(t))
    {

    }
    ~ThreadJoiner()
    {
        if (t.joinable())
            t.join();
    }
private:
    std::thread t;
};


class Range
{
public:

    explicit Range(double begin,double end)
    :begin(begin),end(end)
    {

    }

    double getBegin()const
    {
        return begin;
    }

    double getEnd()const
    {
        return end;
    }

    double getDistance()const
    {
        return end-begin;
    }

private:

double begin;
double end;

};

class Function
{
public:


    Function()=default;
    explicit Function(std::function<double(double)>f
            ,unsigned monteCarloPrecision=defaultPrecisionOfComputing
            , unsigned rectangularIntegrationPrecision=defaultPrecisionOfComputing
            , unsigned trapezoidalIntegrationPrecision=defaultPrecisionOfComputing
            , unsigned simpsonIntegrationPrecision=defaultPrecisionOfComputing)
    :f(f)
            ,monteCarloPrecision(monteCarloPrecision)
    ,rectangularIntegrationPrecision(rectangularIntegrationPrecision)
    ,trapezoidalIntegrationPrecision(trapezoidalIntegrationPrecision)
    ,simpsonIntegrationPrecision(simpsonIntegrationPrecision)
    {

    }

    Range valueRange(Range xRange, unsigned lookingPrecision=defaultPrecisionOfComputing)
    {
        double max=0.0;
        double min=0.0;
        const double singleStep=xRange.getDistance()/static_cast<double>(defaultPrecisionOfComputing);

       for(double funStep=xRange.getBegin();funStep<xRange.getEnd();funStep+=singleStep)
       {
            if(f(funStep)>max)
              max=f(funStep);
            if(f(funStep)<min)
                min=f(funStep);
       }
        return Range(min,max);
    }

    double monteCarlo(Range xRange,Range yrange)
    {
        auto startMonte=std::chrono::system_clock::now();

        std::default_random_engine engine;

        Range yPosRange=Range(0,yrange.getEnd());
        Range yNegRange=Range(yrange.getBegin(),0);

        const double rectangleArea=xRange.getDistance()*yrange.getDistance();
       const double rectanglePosArea=xRange.getDistance()*yPosRange.getDistance();
       const double rectangleNegArea=xRange.getDistance()*yNegRange.getDistance();

        std::uniform_real_distribution<double>xdistribution(xRange.getBegin(),xRange.getEnd());
        std::uniform_real_distribution<double>ydistribution(yrange.getBegin(),yrange.getEnd());

        std::vector<std::pair<double,double>>points(monteCarloPrecision);
        unsigned hittedPointsCounter;
        unsigned negHittedPointsCounter=0;
        unsigned posHittedPointsCounter=0;

      for(auto &vpi:points)
      {
          vpi.first = xdistribution(engine);
          vpi.second = ydistribution(engine);

          if (f(vpi.first)>0? (vpi.second <= f(vpi.first) and vpi.second) : (vpi.second >= f(vpi.first) and
                                                                            !vpi.second))
          {
              vpi.second>0 ? posHittedPointsCounter++ : negHittedPointsCounter++;
          }
      }
        auto endMonte=std::chrono::system_clock::now();

        auto durationMonte=std::chrono::duration<double>(endMonte-startMonte);

        monteCarloDuration=durationMonte;

        hittedPointsCounter=posHittedPointsCounter+negHittedPointsCounter;



        double multiFactor=static_cast<double>(posHittedPointsCounter)/static_cast<double>(hittedPointsCounter)
                           -static_cast<double>(negHittedPointsCounter)/static_cast<double>(hittedPointsCounter);
        std::cout<<"POS"<<posHittedPointsCounter<<"NEG"<<negHittedPointsCounter;

       return multiFactor*static_cast<double>(hittedPointsCounter)*rectangleArea/static_cast<double>(monteCarloPrecision);
    }

    double simpsonIntegration(Range xRange);


    double rectangularIntegration(Range xRange)
    {
        auto startRect=std::chrono::system_clock::now();
        double areaResult=0.0;
        const double rectangleBase=xRange.getDistance()/ static_cast<double>(rectangularIntegrationPrecision);
       for(double xfunStep=xRange.getBegin();xfunStep<=xRange.getEnd();xfunStep+=rectangleBase)
       {
          areaResult+=rectangleBase*f(xfunStep);
       }

        auto endRect=std::chrono::system_clock::now();
        auto durationRect=std::chrono::duration<double>(endRect-startRect);

        rectangularIntegrationDuration=durationRect;

        return areaResult;
    }

    double trapezoidalIntegration(Range xRange)
    {
        auto startTrapezoidal=std::chrono::system_clock::now();
      double areaResult=0.0;
      const double trapezoidBase=xRange.getDistance()/static_cast<double>(trapezoidalIntegrationPrecision);

       for(double xfunStep=xRange.getBegin();xfunStep<=xRange.getEnd();xfunStep+=trapezoidBase)
       {
          areaResult+=trapezoidBase*(f(xfunStep)+f(xfunStep+trapezoidBase))/2;
       }

        auto endTrapezoidal=std::chrono::system_clock::now();
        auto durationTrapezoidal=std::chrono::duration<double>(endTrapezoidal-startTrapezoidal);

        trapezoidalIntegrationDuration=durationTrapezoidal;

        return areaResult;
    }

    double operator()(double x)
    {
     return f(x);
    }

    auto getMonteCarloDuration()const
    {
        return monteCarloDuration;
    }

    auto getRectangularDuration()const
    {
        return rectangularIntegrationDuration;
    }

    auto getTrapezoidalDuration()const
    {
        return trapezoidalIntegrationDuration;
    }

    auto getSimpsonDuration()const
    {
        return simpsonIntegrationDuration;
    }



private:

    std::function<double(double)>f;
    unsigned monteCarloPrecision;
    unsigned rectangularIntegrationPrecision;
    unsigned trapezoidalIntegrationPrecision;
    unsigned simpsonIntegrationPrecision;
    static unsigned defaultPrecisionOfComputing;

    std::chrono::duration<double>monteCarloDuration;
    std::chrono::duration<double>rectangularIntegrationDuration;
    std::chrono::duration<double>trapezoidalIntegrationDuration;
    std::chrono::duration<double>simpsonIntegrationDuration;

};

class SquareFunction:public Function
{
public:


    explicit SquareFunction(double a,double b,double c)
    :a(a),b(b),c(c)
    {

    }

    explicit SquareFunction(Range approxRange,std::function<double(double)>f)
    :approxArea(countDefiniteIntegral(approxRange)),f(f)
    {

    }

    double countDefiniteIntegral(Range integralRange) {
        double a = std::get<0>(countCoefficientsForRange(integralRange));
        double b = std::get<1>(countCoefficientsForRange(integralRange));
        double c = std::get<2>(countCoefficientsForRange(integralRange));


        return (1. / 3) * a * pow(integralRange.getEnd(), 3) + (1. / 2) * b * pow(integralRange.getEnd(), 2) +
               c * integralRange.getEnd()
               - (1. / 3) * a * pow(integralRange.getBegin(), 3) - (1. / 2) * b * pow(integralRange.getBegin(), 2) -
               c * integralRange.getBegin();
    }
    std::tuple<double,double, double> countCoefficientsForRange(Range partRange)
    {
        double mid=(partRange.getBegin()+partRange.getEnd())/2;
        double beg=partRange.getBegin();
        double end=partRange.getEnd();

        double fmid=f(mid);
        double fbeg=f(beg);
        double fend=f(end);

        const double denominator=beg*beg-2*beg*end+end*end;


        c=(beg*beg*fend+beg*end*fend-4*beg*end*fmid+beg*end*fbeg+end*end*fbeg)/denominator;
        b=(-1)*(beg*fbeg-4*beg*fmid+3*beg*fend+3*end*fbeg-4*end*fmid+end*fend)/denominator;
        a=2*(fbeg-2*fmid+fend)/denominator;

      /* std::cout<<mid<<" "<<beg<<" "<<end<<" "<<fmid<<" "<<fbeg<<" "<<fend<<std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
       std::cout<<denominator<<std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
       std::cout<<a<<" "<<b<<" "<<c<<std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));*/

        return std::make_tuple(a,b,c);
    }

    std::pair<double,double>getTop()
    {
        top=std::make_pair(-b/(2*a),f(-b/(2*a)));
        return top;
    }

    double getApproxArea()const
    {
        return approxArea;
    }

private:

    std::function<double(double)>f;
    double a;
    double b;
    double c;
    std::pair<double,double>top;
    double xTop;
    double yTop;
    double delta;
    double approxArea;
};

double fsquare(double x)
{
    return pow(x,2)+2*x+1;
}

double ftrigonometry(double x)
{
    return sin(x)+2*cos(x)+atan(x);
}

double flogarithmic(double x)
{
    return log(x+5);
}

double fexponential(double x)
{
    return exp(2*x);
}

double Function::simpsonIntegration(Range xRange)
{
    auto startSimpson=std::chrono::system_clock::now();
    double areaResult=0.0;
    const double figuresBase=xRange.getDistance()/ static_cast<double>(simpsonIntegrationPrecision);



    for(double xfunStep=xRange.getBegin();xfunStep<=xRange.getEnd();xfunStep+=figuresBase)
    {
        SquareFunction approxSquareFun(Range(xfunStep,xfunStep+figuresBase),f);
        areaResult+=approxSquareFun.getApproxArea();
      //  std::cout<<" "<<areaResult;
      //  std::this_thread::sleep_for(std::chrono::milliseconds(2));

    }
    auto endSimpson=std::chrono::system_clock::now();
    auto durationSimpson=std::chrono::duration<double>(endSimpson-startSimpson);

    Function::simpsonIntegrationDuration=durationSimpson;

    return areaResult;
}


int main()
{

Function squareFun{Function(static_cast<std::function<double(double)>>(fsquare))};
Function trigoFun{Function(static_cast<std::function<double(double)>>(ftrigonometry))};
Function logFun{Function(static_cast<std::function<double(double)>>(flogarithmic))};
Function expFun{Function(static_cast<std::function<double(double)>>(fexponential))};

//Range userDefinedRange{Range(5.0,11.0)};

    double userDefinedBeg;
    double userDefinedEnd;


    std::cout<<"input range begin x coordinate: ";
    std::cin>>userDefinedBeg;
    std::cout<<std::endl;

    std::cout<<"input range end x coordinate: ";
    std::cin>>userDefinedEnd;
    std::cout<<std::endl;


    Range userDefinedRange{Range(userDefinedBeg,userDefinedEnd)};

    std::vector<std::tuple<Function,std::string,double,double,double,double>>functionsWithAreasMulti;
    std::vector<std::tuple<Function,std::string,double,double,double,double>>functionsWithAreasSingle;

    using Type=decltype(functionsWithAreasMulti);

    auto squareCompute=[&](Type &con)
    {
                con.emplace_back(std::make_tuple(
                squareFun
                ,GET_VARIABLE_NAME(squareFun)
                ,squareFun.monteCarlo(userDefinedRange,squareFun.valueRange(userDefinedRange))
                ,squareFun.rectangularIntegration(userDefinedRange)
                ,squareFun.trapezoidalIntegration(userDefinedRange)
                ,squareFun.simpsonIntegration(userDefinedRange)));
    };

    auto trigCompute=[&](Type &con)
    {
                con.emplace_back(std::make_tuple(
                trigoFun, GET_VARIABLE_NAME(trigoFun),
                trigoFun.monteCarlo(userDefinedRange, trigoFun.valueRange(userDefinedRange)),
                trigoFun.rectangularIntegration(userDefinedRange), trigoFun.trapezoidalIntegration(userDefinedRange),
                trigoFun.simpsonIntegration(userDefinedRange)));
    };

    auto logCompute=[&](Type &con)
    {
                con.emplace_back(std::make_tuple(
                logFun, GET_VARIABLE_NAME(logFun), logFun.monteCarlo(userDefinedRange, logFun.valueRange(userDefinedRange)),
                logFun.rectangularIntegration(userDefinedRange), logFun.trapezoidalIntegration(userDefinedRange),
                logFun.simpsonIntegration(userDefinedRange)));
    };

    auto expCompute=[&](Type &con)
    {
        con.emplace_back(std::make_tuple(
                expFun, GET_VARIABLE_NAME(expFun), expFun.monteCarlo(userDefinedRange, expFun.valueRange(userDefinedRange)),
                expFun.rectangularIntegration(userDefinedRange), expFun.trapezoidalIntegration(userDefinedRange),
                expFun.simpsonIntegration(userDefinedRange)));
    };

    std::thread squareThread{squareCompute,std::ref(functionsWithAreasMulti)};
    std::thread trigoThread{trigCompute,std::ref(functionsWithAreasMulti)};
    std::thread logThread{logCompute,std::ref(functionsWithAreasMulti)};
    std::thread expThread{expCompute,std::ref(functionsWithAreasMulti)};

    auto startMulti=std::chrono::system_clock::now();

    squareThread.join();
    trigoThread.join();
    logThread.join();
    expThread.join();


    auto endMulti=std::chrono::system_clock::now();

    auto durationMulti=std::chrono::duration<double>(endMulti-startMulti);



    for(const auto &vti:functionsWithAreasMulti)
{
    std::cout<<std::endl;
    std::cout<<"MULTITHREADING WAY: ";
    std::cout<<std::endl;
    std::cout<<std::get<1>(vti)<<std::endl<<std::endl<<"Monte carlo integration: "
             <<std::get<2>(vti)
             <<std::endl<<"computations in this method took: "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getMonteCarloDuration()).count() <<" miliseconds"
             <<std::endl<<std::endl<<"Rectangle integration: "<<std::get<3>(vti)
             <<std::endl<<"computations in this method took: "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getRectangularDuration()).count() <<" miliseconds"
             <<std::endl<<std::endl<<"Trapezoid integration: " <<std::get<4>(vti)
             <<std::endl<<"computations in this method took: "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getTrapezoidalDuration()).count() <<" miliseconds"
             <<std::endl<<std::endl<<"Simpson integration: "<<std::get<5>(vti)
             <<std::endl<<"computations in this method took:"
             <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getSimpsonDuration()).count() <<" miliseconds"
             <<std::endl <<std::endl <<"in total computations of this function definite integral took: " <<std::chrono::duration_cast<std::chrono::milliseconds>(durationMulti).count()
             <<" miliseconds"<<std::endl;
}



    auto startSingle=std::chrono::system_clock::now();

    squareCompute(functionsWithAreasSingle);
    trigCompute(functionsWithAreasSingle);
    logCompute(functionsWithAreasSingle);

    auto endSingle=std::chrono::system_clock::now();
    auto durationSingle=std::chrono::duration<double>(endSingle-startSingle);



    for(const auto &vti:functionsWithAreasSingle)
    {
        std::cout<<std::endl;
        std::cout<<"SINGLE-THREAD WAY: ";
        std::cout<<std::endl;
        std::cout<<std::get<1>(vti)<<std::endl<<std::endl<<"Monte carlo integration: "
                 <<std::get<2>(vti)
                 <<std::endl<<"computations in this method took: "
                 <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getMonteCarloDuration()).count() <<" miliseconds"
                 <<std::endl<<std::endl<<"Rectangle integration: "<<std::get<3>(vti)
                 <<std::endl<<"computations in this method took: "
                 <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getRectangularDuration()).count() <<" miliseconds"
                 <<std::endl<<std::endl<<"Trapezoid integration: " <<std::get<4>(vti)
                 <<std::endl<<"computations in this method took: "
                 <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getTrapezoidalDuration()).count() <<" miliseconds"
                 <<std::endl<<std::endl<<"Simpson integration: "<<std::get<5>(vti)
                 <<std::endl<<"computations in this method took:"
                 <<std::chrono::duration_cast<std::chrono::milliseconds>(std::get<0>(vti).getSimpsonDuration()).count() <<" miliseconds"
                 <<std::endl <<std::endl<<"in total computations of this function definite integral took: " <<std::chrono::duration_cast<std::chrono::milliseconds>(durationSingle).count()
                 <<" miliseconds"<<std::endl;

    }

    return 0;
}

unsigned Function::defaultPrecisionOfComputing=1000000;
