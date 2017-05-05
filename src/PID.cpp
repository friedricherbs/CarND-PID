#include "PID.h"
#include <fstream>      // std::ofstream
#include <cmath>
#include <assert.h>     /* assert */

using namespace std;

/*
* TODO: Complete the PID class.
*/
#define TOL               0.2
#define INIT_PHASE_LENGTH 100
#define ACCU_PHASE_LENGTH 4000

PID::PID() {}

PID::~PID() {}

void PID::Init( const double Kp,  const double Ki,  const double Kd, 
                const double dKp, const double dKi, const double dKd,
				const unsigned int twiddleUpdateParam, 
				const double       bestErr, 
				const TwiddleMode  mode                                 ) 
{
  m_Kp                 = Kp;
  m_Ki                 = Ki;
  m_Kd                 = Kd;
                      
  m_p_error            = 0.0;
  m_i_error            = 0.0;
  m_d_error            = 0.0;
                      
  m_err                = 0.0;
  m_dp_p               = dKp;
  m_dp_i               = dKi;
  m_dp_d               = dKd;
  
  m_twiddleUpdateParam = twiddleUpdateParam; 
  m_bestError          = bestErr;
  
  m_twiddleMode        = mode;

  m_cycleCounter       = 0U;
}

void PID::UpdateError(const double cte) 
{
  m_d_error  = cte - m_p_error;
  m_p_error  = cte; 
  m_i_error += cte;

  ++m_cycleCounter;
  
  if (m_cycleCounter > INIT_PHASE_LENGTH)
  {
	  m_err += std::pow(cte, 2.0);
  }
}

double PID::TotalError() const
{
  return -m_Kd*m_d_error - m_Kp*m_p_error - m_Ki*m_i_error;
}


bool PID::Twiddle() 
{
	if (m_cycleCounter < ACCU_PHASE_LENGTH)
	{
	  return false;
	}
	
    const double delta = m_dp_d + m_dp_i + m_dp_p;
	if (delta > TOL)
	{
		switch (m_twiddleMode)
		{
		case TWIDDLE_INIT:
		  {
			TwiddleReInit(m_err, 0.0, 1.0, 1.0, m_twiddleUpdateParam, TWIDDLE_TRIAL);
		  }
		  break;
		case TWIDDLE_TRIAL:
		  {
		    if (m_err < m_bestError)
			{ 
		      const unsigned int nextUpdateParam = (m_twiddleUpdateParam+1)%3;
		      TwiddleReInit(m_err, 1.0, 1.1, 0.0, nextUpdateParam, TWIDDLE_TRIAL);
			}
			else
			{
			  TwiddleReInit(m_bestError, 0.0, 1.0, -2.0, m_twiddleUpdateParam, TWIDDLE_CHANGE_SIGN);
			}
		  }
		  break;
		case TWIDDLE_CHANGE_SIGN:
		  {
			const unsigned int nextUpdateParam = (m_twiddleUpdateParam+1)%3;
		    if (m_err < m_bestError)
			{ 
		      TwiddleReInit(m_err, 1.0, 1.1, 0.0, nextUpdateParam, TWIDDLE_TRIAL);
			}
			else
			{
		      TwiddleReInit(m_bestError, 1.0, 0.9, 1.0, nextUpdateParam, TWIDDLE_TRIAL);
			}
		  }
		  break;
		default:
		  {
		    assert(0 && "Unknown twiddle mode");
		  }
		  break;			
		}
		
		return true;
	}
	return false;
}

void PID::TwiddleReInit(const double best_err,                const double pChangeValFactorNext,  const double dpIncreaseFactor, 
                        const double pChangeValFactorCurrent, const unsigned int nextUpdateParam, const TwiddleMode nextMode)
{
  // Get current twiddle parameters
  double p[3]  = {m_Kp,   m_Ki,   m_Kd};
  double dp[3] = {m_dp_p, m_dp_i, m_dp_d};
  
  // Dump last parameters
  std::ofstream ofs;
  ofs.open("params.txt", std::ofstream::out | std::ofstream::app);
  ofs << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << dp[0] << "\t" << dp[1] << "\t" << dp[2] << "\t" << m_err << "\t" << m_twiddleUpdateParam << "\t" << m_twiddleMode << "\t" << nextUpdateParam << "\n";
  ofs.close();
  
  // Change parameter for current coordinate
  p[m_twiddleUpdateParam]  += (dp[m_twiddleUpdateParam]*pChangeValFactorCurrent);

  // Change step size for this coordinate
  dp[m_twiddleUpdateParam] *= dpIncreaseFactor;  
  
  // Change parameter for next run
  p[nextUpdateParam]       += (dp[nextUpdateParam]*pChangeValFactorNext);
  
  // Set new filter parameters
  Init(p[0], p[1], p[2], dp[0], dp[1], dp[2], nextUpdateParam, best_err, nextMode);
}

