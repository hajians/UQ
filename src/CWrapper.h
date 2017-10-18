#ifndef __CWRAPPER_H
#define __CWRAPPER_H

#ifdef __cplusplus
  extern "C"
  {
#endif

    typedef struct SemilinearSystem SemilinearSystem;

    SemilinearSystem* CSemiLinSystem(double c_sound, double t_final,
				     double x_l, double x_r,
				     double dx, int lambda_len);

    void CRun(SemilinearSystem* pipe, double coef[]);
    void CWrite2File(SemilinearSystem* pipe, char* filename, bool append);
    void CInfo(SemilinearSystem* pipe);

    int CCurrentTimeIndex(SemilinearSystem* pipe);

    double* CBoundaryValueP_Left(SemilinearSystem* pipe);
    double* CBoundaryValueP_Right(SemilinearSystem* pipe);
      
#ifdef __cplusplus
  }
#endif

#endif
