#include "SemilinearSystem.h"
#include "CWrapper.h"

extern "C" {

  SemilinearSystem* CSemiLinSystem(double c_sound, double t_final,
				   double x_l, double x_r,
				   double dx, int lambda_len,
				   double eps)
  {
    return new SemilinearSystem(c_sound, t_final, x_l, x_r, dx, lambda_len, eps);
  }

  void CRun(SemilinearSystem* pipe, double coef[], bool write2file_bool)
  {
    pipe->Run(coef, write2file_bool);
  }
  void CInfo(SemilinearSystem* pipe)
  {
    pipe->info();
  }

  void CWrite2File(SemilinearSystem* pipe, char* filename, bool append)
  {
    pipe->Write2File(filename, append);
  }

  int CCurrentTimeIndex(SemilinearSystem* pipe)
  {
    return pipe->CurrentTimeIndex();
  }

  double* CBoundaryValueP_Left(SemilinearSystem* pipe)
  {
    return pipe->BoundaryValueP_Left();
  }

  double* CBoundaryValueP_Right(SemilinearSystem* pipe)
  {
    return pipe->BoundaryValueP_Right();
  }

  double* CTimeSlices(SemilinearSystem* pipe)
  {
    return pipe->TimeSlices();
  }

}
