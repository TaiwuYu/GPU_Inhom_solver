 #ifndef CONSTANTS_H
 #define CONSTANTS_H
 // Constants for the HEA model
 
//#define BIN_FILE

//#define FROM_INPUT_FILE
  //  	#define RESTART_FROM_END
		#define RESTART_TIMESTEP 0

#define SAVE_CONFIG 1000 

//#define NOISE_ON

//#define ALLOW_APB

 
 //Constants in the free energy
 #define A1 (0.5)
 #define A2 (1.0)
 #define B1 (0.1)
 #define W	(0.0)
 #define V  (1)
 //Constants of the system size, discretization
 #define L (128)
   #define M (128)
      #define N (128)
      #define LMN (L*M*N)
#define HLMN (L*M*(N/2+1))
 #define dx (1.0)
 #define dy (1.0)
 #define dz (1.0)
 #define dt_value (0.02)

//Kinetic parameters and gradient coefficient
 #define Mc 	(1.0)
 #define Meta 	(5.0)
 #define kc    (0.6)
 #define keta  (6.0)

//Noise parameters
	#define Namp (2.0)
	#define Namp2 (0.0)
	#define NRAND (1000)

//#define DEBUG_MODE_ON

#define MAXITER (100) // maximum iterations for the inhomogeneous solver
//#define MISFIT_COUPLED_WITH_CONC
//#define MISFIT_COUPLED_WITH_ETA
#endif

  
 
 
 
