#include "include/cxxopts.hpp"

cxxopts::ParseResult parse(int argc, char *argv[]) {
  try {
    cxxopts::Options options(argv[0], " -testing cxxopts");
    options.positional_help("[optional args]").show_positional_help();

    options.allow_unrecognised_options().add_options()(
  "l, length", "length of one side of the domain.",
  cxxopts::value<double>(),
  "l")("c, coeff", "coefficient of the initialization function.",
       cxxopts::value<double>(), "c")(
  "n, nu", "nu value (viscosity) to use", cxxopts::value<double>(), "nu")(
  "e, eta", "Eta value (a parameter in the initial condition)",
  cxxopts::value<double>(), "eta")(
  "g, gradVratio_bound", "Gradient of Vratio upper bound. One of the criteria for smoothing.",
  cxxopts::value<double>(), "gradVratio_bound")(
  "b, L2ratio_bound", "Upper bound of ratio of L2. One of the criteria for smoothing. ",
  cxxopts::value<double>(), "L2ratio_bound")("p, print", "Print every p-th step.",
           cxxopts::value<int>(), "print")(
  "f, frequency", "How often apply coarse-fine.", cxxopts::value<int>(),
  "freq")(
"o, cfonoff", "Turn coarse fine on (1) or off(0).", cxxopts::value<int>(),
"cfonoff")("s, silent", "do not print calculations on screen");

    // options.parse_positional({"input", "output", "positional"});
    auto result = options.parse(argc, argv);

    // std::cout << "Arguments remain = " << argc << std::endl;
    return result;
  } catch (const cxxopts::OptionException &e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

