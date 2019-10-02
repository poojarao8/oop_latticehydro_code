#include "main.h"

using namespace std;

int main(int argc, char *argv[]) {
  int outdatasize;
  // values to report:
  // L,W,H,mu,eta,nu,h,t,VnormL2,VnormLinf,gradVnormL2,pct1,pct2,divVLinf
  double *outdata = new double[10001 * 16];
  double VnormL2, maximumV, pct1, pct2, divVLinf, gradnorm, gradsymm, gradskew, gradnorm0;
  double maximumVC, gradnormC, gradsymmC, gradskewC;
  double *V = new double[ARR_SIZE];		// Velocity on fine grid
  double *Vc = new double[COARSE_ARR_SIZE];     // coarse grid velocity
  double *tempC = new double[COARSE_ARR_SIZE];  // coarse grid temp arr.
  double *Laplace_V = new double[ARR_SIZE];
  double *Vfvf_res = new double[ARR_SIZE];
  double *temp = new double[ARR_SIZE];
  double *R = new double[ARR_SIZE];
  double *Total = new double[ARR_SIZE];
  double *scal = new double[L * W * H];
  double *pressure_old = new double[L * W * H];
  double *solution_old = new double[L * W * H];
  double *gradnormdomain = new double[L * W * H];
  double *gradnormdomainC = new double[LL*WW*HH];
  double threshold1, threshold2;
  double dt, t;
  double nu = 0.001;
  double eta = 100.0;
  double h = 0.0;
  double phys_length = 6.283185307179586;
  double coeff = MU;
  double L2ratio_bound =1.4;
  double gradVratio_bound=1.4;
  double V0L2;
  double gradV0ratio, gradVratio;
  double smoothing_reason=0.0; //Reason for smoothing applied. 0=not applied, 1=L2, 2=gradVratio, 3=periodic
  int pr = 1;
  int freq = 1000000; //For periodic application of smoothing. Set freq=1000000 to turn this off.
  int basetime = 0;
  auto result = parse(argc, argv);
  auto arguments = result.arguments();
  bool silent = false;
  int cfonoff = 1;
  int CF_FLAG=0;
  int smoothing_count = 0;

  if (result.count("c")) coeff = result["c"].as<double>();
  if (result.count("l")) phys_length = result["l"].as<double>();
  if (result.count("nu")) nu = result["nu"].as<double>();
  if (result.count("eta")) eta = result["eta"].as<double>();
  if (result.count("silent")) silent = true;
  if (result.count("print")) pr = result["print"].as<int>();
  if (result.count("freq")) freq = result["freq"].as<int>();
  if (result.count("cfonoff")) cfonoff = result["cfonoff"].as<int>();
  if (result.count("L2ratio_bound")) L2ratio_bound = result["L2ratio_bound"].as<double>();
  if (result.count("gradVratio_bound")) gradVratio_bound = result["gradVratio_bound"].as<double>();

  if (!silent) {
    std::cout << "Using values : phys_length " << phys_length << ", nu " << nu
	      << ", eta " << eta << ", smoothing frequency:" << freq
	      << std::endl;
  }

  h = phys_length / L;

  // Testing
  char plane[] = "xyz";
  for (int i = 0; i < COARSE_ARR_SIZE; i++) Vc[i] = 0.0;
  for (int i = 0; i < (L * W * H); i++) scal[i] = 0.0;
  for (int i = 0; i < (L * W * H); i++) pressure_old[i] = 0.0;
  for (int i = 0; i < (L * W * H); i++) solution_old[i] = 0.0;

  outdatasize = 0;
  t = 0.0;
  dt = 0.01;

  if (nu!=0) {
    if (dt>0.5*h*h/nu) {
       cout<<RED<<"dt too large compared to h^2/nu. Adjusting."<<RESET<<endl;
       dt = 0.5*h*h/nu;
    }
  }

  cout << "dt: " << dt << endl;
  for (int i = 0; i < COARSE_ARR_SIZE; i++) tempC[i] = 0.0;

  // initialize a divergence free vector field.
  initialize_var(tempC, h, eta, coeff);
  bdC(tempC, Vc);  // This makes Vc a divergence free vector field since
		   // div(curl (vector)) = 0.

  coarse2fine(Vc, V, WW, LL, HH);  // make V the fine version of Vc

  cout << "L_inf norm of div V is: " << maxnormscal(scal, SCL_ARR_SIZE) << endl;
  cout << "Lattice Size: (" << L << "," << W << ","
       << "," << H << ")" << endl;
  cout << "h: " << h << endl;
  cout << "Viscosity: " << nu << endl;

  V0L2 = norm(V, ARR_SIZE);
  threshold1 = 0.01 * maxnorm(V, ARR_SIZE);
  //gradV0ratio = grad_ratio_norm(V);
  gradient_norm(V, L, W, H, &gradnorm0, &gradsymm, &gradskew);

  //calculate the L2 norm of the grad
  gradient_norm_domain(V, L, W, H, gradnormdomain);

  //calculate the L2 norm of the grad on coarse grid
  gradient_norm_domain(Vc, LL, WW, HH, gradnormdomainC);

  /* adding code for writing visit files at initial time before RK4*/

  char outfile[100];
  output_file(0, nu, outfile);
  //visit_plot(outfile, L, W, H, pressure_old, V);
  visit_plot(outfile, LL, WW, HH, gradnormdomainC, Vc, t, smoothing_count);

  for (int rep = 0; rep < 6001; rep++) {
    maximumV = maxnorm(V, ARR_SIZE);
    if (maximumV != maximumV)
    {  // if maximumV is nan.
      cout << "Blowup encountered." << endl;
      break;
    }

    if (maximumV < 1e-6)
    {
      cout << "Decays to zero." << endl;
      break;
    }
    if (dt>0.5*h/maximumV) {
       cout<<RED<<"dt too large compared to h/umax. Adjusting."<<RESET<<endl;
       dt = 0.5*h/maximumV;
    }

    VnormL2 = norm(V, ARR_SIZE);
    gradient_norm(V, L, W, H, &gradnorm, &gradsymm, &gradskew);

    threshold2 = 0.25 * maximumV;
    pct_vectors(V, ARR_SIZE, L, W, H, threshold1, threshold2, &pct1, &pct2);
    bd10(V, scal);
    divVLinf = maxnormscal(scal, SCL_ARR_SIZE);
    //gradVratio = grad_ratio_norm(V);

    // check cf_flag here
    CF_FLAG=0;
    smoothing_reason=0.0;

    if (VnormL2/V0L2 > L2ratio_bound) {CF_FLAG = 1; smoothing_reason=1.0;
    cout<<RED << "L2 norm of velocity high. Smoothing." <<RESET <<endl;
     }
    if ( gradnorm/gradnorm0 > gradVratio_bound ) {CF_FLAG = 1; smoothing_reason=2.0;
      cout<<RED << "Gradient of V high. Smoothing." <<RESET <<endl;
    }

    if((rep+1)%freq==0) {CF_FLAG=1; smoothing_reason=3.0; }
    //

    if ((cfonoff==1)&&(CF_FLAG==1)) {
      cout << "Smoothing done." <<endl;
      fine2coarse(V, Vc, L, W, H);
      coarse2fine(Vc, V, LL, WW, HH);
      V0L2=norm(V, ARR_SIZE);
      gradient_norm(V, L, W, H, &gradnorm0, &gradsymm, &gradskew);
      //gradV0ratio=grad_ratio_norm(V);
      smoothing_count += 1;
    }

    if ((rep+1) % pr == 0) {

      if (!silent) {
	      cout << "h,nu: " << h << " , " << nu << endl;
	      cout << "time: " << t << endl;
	      cout << "step: " << rep << endl;
        cout << "smoothing count " << smoothing_count << endl;
	      cout << "velocity L-inf norm: " << maximumV << endl;
	      cout << "velocity L2 norm: " << VnormL2 << endl;
	      cout << "grad(V) L2 norm: " << gradnorm << endl;
	      cout << "L2 norm of symmetric part of grad(V): " << gradsymm << endl;
	      cout << "L2 norm of skew-symmetric part of grad(V): " << gradskew
	           << endl;
	      cout << "Percent of velocities with norm greater than " << threshold1
	           << " :" << pct1 << endl;
	      cout << "Percent of velocities with norm greater than " << threshold2
	           << " :" << pct2 << endl;
        //cout << "Gradient Vratio : " << gradVratio << endl;
	      cout << "divergence L-inf norm: " << divVLinf << endl;
	      cout << endl;
      }

      // values to report:
      // L,W,H,mu,eta,nu,h,t,VNormL2,VnormLinf,gradnorm,gradsymm,gradskew,pct1,pct2,divVLinf
      outdata[16 * outdatasize + 0] = L;
      outdata[16 * outdatasize + 1] = W;
      outdata[16 * outdatasize + 2] = H;
      outdata[16 * outdatasize + 3] = MU;
      outdata[16 * outdatasize + 4] = eta;
      outdata[16 * outdatasize + 5] = nu;
      outdata[16 * outdatasize + 6] = h;
      outdata[16 * outdatasize + 7] = t;
      outdata[16 * outdatasize + 8] = VnormL2;
      outdata[16 * outdatasize + 9] = maximumV;
      outdata[16 * outdatasize + 10] = gradnorm;
      outdata[16 * outdatasize + 11] = gradsymm;
      outdata[16 * outdatasize + 12] = gradskew;
      outdata[16 * outdatasize + 13] = pct1;
      outdata[16 * outdatasize + 14] = pct2;
      outdata[16 * outdatasize + 15] = smoothing_reason;
      outdatasize++;
    }

    // Runge Kutta
    nav_stoke(h, nu, V, Vfvf_res, R, scal, pressure_old,
	      solution_old);  // scal is pressure.
    for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt * 0.5;
    for (int i = 0; i < ARR_SIZE; i++)
      Total[i] = R[i];  // separate loop to make compiler to use SIMD
    nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
    for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt * 0.5;
    for (int i = 0; i < ARR_SIZE; i++) Total[i] += 2 * R[i];
    nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
    for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt;
    for (int i = 0; i < ARR_SIZE; i++) Total[i] += 2 * R[i];
    nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
    for (int i = 0; i < ARR_SIZE; i++) Total[i] += R[i];
    for (int i = 0; i < ARR_SIZE; i++) V[i] += Total[i] * (dt / 6);
    t += dt;


    // print out maxnormC and gradient every 50th step
    if ((rep+1) % 10 == 0) {
      fine2coarse(V, Vc, L, W, H);
      maximumVC = maxnorm(Vc, COARSE_SCL_ARR_SIZE);
      gradient_norm(Vc, L, W, H, &gradnormC, &gradsymmC, &gradskewC);

      cout << "Coarse scale velocity L-inf norm: " << maximumVC << endl;

      cout << "Coarse scale grad(V) L2 norm: " << gradnormC << endl;
      cout << "Coarse scale L2 norm of symmetric part of grad(V): " << gradsymmC
	         << endl;
      cout << "Coarse scale L2 norm of skew-symmetric part of grad(V): "
	         << gradskewC << endl
	   << endl;
     write2file(outdata,outdatasize,nu,eta,cfonoff,L2ratio_bound,gradVratio_bound);
    }

    /* adding code for writing visit files */
    if ((rep+1)%20==0)
    {
      char outfile[100];
      output_file(rep+1, nu, outfile);
      fine2coarse(V, Vc, L, W, H);
      gradient_norm_domain(Vc, LL, WW, HH, gradnormdomainC);
      //visit_plot(outfile, L, W, H, pressure_old, V);
      visit_plot(outfile, LL, WW, HH, gradnormdomainC, Vc, t, smoothing_count);
    }
  }

   write2file(outdata,outdatasize,nu,eta,cfonoff,L2ratio_bound,gradVratio_bound);

  return 1;
}
