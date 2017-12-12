#include "SemilinearSystem.h"
#include "CWrapper.h"

extern "C" {

  /**
   * Instantiate an object of SemilinearSystem class.
   */
  SemilinearSystem* CSemiLinSystem(double c_sound, double t_final,
				   double x_l, double x_r,
				   double dx, int lambda_len,
				   double eps)
  {
    SemilinearSystem* s = new SemilinearSystem(c_sound, t_final, x_l, x_r, dx, lambda_len, eps);
    cout << "Address of the pipe in the memory from cpp: "<< s << endl;
    return s;
  }

  /**
   * Computes the solution of the semilinear system.
   */
  void CRun(SemilinearSystem* pipe, double coef[], bool write2file_bool, bool progress_bool)
  {
    pipe->Run(coef, write2file_bool, progress_bool);
  }

  /**
   * Get the info of the pipe.
   */
  void CInfo(SemilinearSystem* pipe)
  {
    pipe->info();
  }

  /**
   * Writes the solution computed using the Run command into a file.
   */
  void CWrite2File(SemilinearSystem* pipe, char* filename, bool append)
  {
    pipe->Write2File(filename, append);
  }

  /**
   * Returns the current time index of the semilinear system.
   */
  int CCurrentTimeIndex(SemilinearSystem* pipe)
  {
    return pipe->CurrentTimeIndex();
  }

  /**
   * Returns the value of P on the left boundary.
   */
  double* CBoundaryValueP_Left(SemilinearSystem* pipe)
  {
    return pipe->BoundaryValueP_Left();
  }

  /**
   * Returns the value of P on the right boundary.
   */
  double* CBoundaryValueP_Right(SemilinearSystem* pipe)
  {
    return pipe->BoundaryValueP_Right();
  }

  /**
   * Gives back an array containing the time slices.
   */
  double* CTimeSlices(SemilinearSystem* pipe)
  {
    return pipe->TimeSlices();
  }

  /**
   * Gives back the Friction function averaged over the mesh.
   */
  double* CLambda_Average(SemilinearSystem* pipe, double DA_P_Lambda_Coefficients_GIVEN[])
  {
    return pipe->LambdaAverage(DA_P_Lambda_Coefficients_GIVEN);
  }

  /**
   * Gets the current friction function averaged over the mesh. This
   * is updated after each Run command.
   */
  double* CGetLambda_Average(SemilinearSystem* pipe)
  {
    return pipe->SendLambdaAverage();
  }
  /**
   * Gives back the number of cells.
   */
  int CNumberofCells(SemilinearSystem* pipe)
  {
    return pipe->NumberofCells();
  }
}
