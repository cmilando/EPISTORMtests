import PackageABC from 'packageABC'
import PackageXYZ from 'packageXYZ'

interface IPackage {
  doCalculation(input: object): number;
}

class Calculator {
  // each calculator has a package
  // this is private
  private package:IPackage;
  
  // depending on the package, you have other options
  constructor(options: {
    package: IPackage;
    anythingElseYouNeed?: string;
  }) {
    this.package = options.package;
  }
  
  // this is the common output
  public doCalculation(input: object): number {
    return this.package.doCalculation(input);
  }
  
  // this the input function for package ABC
  public static fromPackageABC(){
    // some consts
    const abc = new PackageABC({
      whatever:123
    })
    // and creates a new object
    return new Calculator({
      package: {
        doCalculation: (i) => abc.runGaussianRegressionWhatever(i)
      }
    })
  }

  // this the input function for package XYZ
  public static fromPackageXYZ(){
    // some consts
    const xyz = new PackageXYZ({
      sigFigs:1337
    })
    return new Calculator({
      package: {
        doCalculation: (i) => xyz.whoEvenWasBernoulli(i)
      }
    })
  }
}

const calculatorAbc = Calculator.fromPackageABC();
const calculatorXyz = Calculator.fromPackageXYZ();
