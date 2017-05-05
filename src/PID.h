#ifndef PID_H
#define PID_H

class PID {
public:
  /*
  * Errors
  */
  double m_p_error;
  double m_i_error;
  double m_d_error;

  /*
  * Coefficients
  */ 
  double m_Kp;
  double m_Ki;
  double m_Kd;

  /*
  * Algo cycle counter.
  */
  unsigned int m_cycleCounter;
  
  /*
  * Accumulated squared cte for Twiddle.
  */
  double m_err;
  
  /*
  * Twiddle Update parameters.
  */
  double m_dp_p;
  double m_dp_i;
  double m_dp_d;
  
  unsigned int m_twiddleUpdateParam;
  
  double m_bestError;
  
  enum TwiddleMode
  {
	  TWIDDLE_INIT,
	  TWIDDLE_TRIAL,
	  TWIDDLE_CHANGE_SIGN
  };
  
  TwiddleMode m_twiddleMode;

  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();

  /*
  * Initialize PID.
  */
  void Init(const double Kp,  const double Ki,  const double Kd, 
            const double dKp, const double dKi, const double dKd,
			const unsigned int twiddleUpdateParam, 
		    const double bestErr, 
			const TwiddleMode mode  );

  /*
  * Update the PID error variables given cross track error.
  */
  void UpdateError(const double cte);

  /*
  * Calculate the total PID error.
  */
  double TotalError() const;
  
  /*
  * Do twiddle parameter update. Return if one update cycle was done or not.
  */
  bool Twiddle();
  
  private:
  
  /*
  * Set twiddle parameters for next cycle.
  */
  void TwiddleReInit( const double       best_err, 
                      const double       pChangeValFactorNext, 
					  const double       dpIncreaseFactor, 
				      const double       pChangeValFactorCurrent, 
					  const unsigned int nextUpdateParam, 
					  const TwiddleMode  nextMode);
};

#endif /* PID_H */
