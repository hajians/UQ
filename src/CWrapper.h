#ifndef __CWRAPPER_H
#define __CWRAPPER_H

#ifdef __cplusplus
  extern "C"
  {
#endif

    typedef struct SemilinearSystem SemilinearSystem;

    SemilinearSystem* CSemiLinSystem(double c_sound, double t_final,
				     double x_l, double x_r,
				     double dx, int lambda_len, double eps);

    void CRun(SemilinearSystem* pipe, double coef[], bool write2file_bool);

    void CWrite2File(SemilinearSystem* pipe, char* filename, bool append);

    void CInfo(SemilinearSystem* pipe);

    int CCurrentTimeIndex(SemilinearSystem* pipe);

    double* CBoundaryValueP_Left(SemilinearSystem* pipe);

    double* CBoundaryValueP_Right(SemilinearSystem* pipe);

    double* CTimeSlices(SemilinearSystem* pipe);

    double* CLambda_Average(SemilinearSystem* pipe, double DA_P_Lambda_Coefficients_GIVEN[]);

    double* CGetLambda_Average(SemilinearSystem* pipe);
    
    int CNumberofCells(SemilinearSystem* pipe);

#ifdef __cplusplus
  }
#endif

#endif
