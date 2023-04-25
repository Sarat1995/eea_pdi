"""model.py"""

import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.integrate import solve_ivp
import bottleneck as bn

# Biexponential fitting function.
def biexponential(x, a, b, c, d):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        return a * np.exp(b * x) + c * np.exp(d * x)

class Model:

    def __init__(self, b, d, n_start,suppress_excimer_annihilation):

        # Rate constants.
        self.k = np.array([b, d])
        
        # Initial exciton densities at each segment.
        self.n_start = n_start

        self.segment_idx = 0

        self.suppress_excimer_annihilation = suppress_excimer_annihilation

        self.fitting_power = 0

    # Evaluate model TRPL data for given parameters.
    def trpl_model_segment(self, t, gamma, update_n_start=True):

        N = self.n_start.shape[1]

        if self.fitting_power == 0:
            fitting_powers = range(N)
        else:
            fitting_powers = [0, self.fitting_power]

        # Array to store model TRPL data.
        model_data = np.zeros((len(fitting_powers), len(t)))

        # For each power:
        for i in range(len(fitting_powers)):

            # Initial populations of two species.
            n0 = self.n_start[self.segment_idx, fitting_powers[i]]# * n0_relative[i]
    
            # dn/dt function.
            def deriv(t, n):

                assert self.k[0]>self.k[1],"Components not ordered correctly with faster component first."
                
                if self.suppress_excimer_annihilation:
                    return -self.k * n - np.array([gamma*n[0]*np.sum(n),0])
                else:
                    return -self.k * n - gamma*n*np.sum(n)
                    
            # Solve model numerically.
            sol = solve_ivp(deriv, (t[0], t[-1]), n0, t_eval=t)
            assert sol.success
            model_data[i] = np.sum(sol.y, axis=0)
    
            if update_n_start:
                self.n_start[self.segment_idx + 1, fitting_powers[i]] = sol.y[:, -1]

        return model_data

    def trpl_model(self, t_segments, gamma):

        num_segments = len(t_segments)

        result = []

        self.segment_idx = 0

        for i in range(num_segments):

            t_segment = t_segments[i]

            segment_result = self.trpl_model_segment(t_segment, gamma[i])

            if i == num_segments - 1:
                result.append(segment_result)
            else:
                result.append(segment_result[:, :-1])

            self.segment_idx += 1

        self.segment_idx = 0

        result = np.concatenate(result, axis=1)

        return result

    # Perform log least-squares fit for gamma.
    def fit_scaled_gamma(self, t_segments, smoothed_relative_segments):

        num_segments = len(t_segments)

        output = []
        error = []

        # Define logarithmic fitting function.
        def log_model(t, gamma):#, n_start=np.array([a, c])):
            return (self.trpl_model_segment(t, gamma)).ravel()

        self.segment_idx = 0

        for i in range(num_segments):

            t_segment = t_segments[i] 
            smoothed_relative_segment = smoothed_relative_segments[i]

            popt, pcov = curve_fit(log_model, t_segment, (smoothed_relative_segment.ravel()), bounds=(0, 10), absolute_sigma=True)

            output.append(popt[0])
            error.append(np.sqrt(pcov[0]))

            self.segment_idx += 1

        self.segment_idx = 0
        #print("gamma = " + str(popt[0] / powers[0]))
    
        output = np.array(output)
        error = np.array(error)

        return output, error

    def maximize_scaled_gamma(self, t, smoothed_relative, margin):

        constraints = []

        fun_above = lambda g:  (1 + margin)*smoothed_relative.ravel() - self.trpl_model([t], [g]).ravel()
        fun_below = lambda g: -(1 - margin)*smoothed_relative.ravel() + self.trpl_model([t], [g]).ravel()
        #fun_above = lambda g: 2*bn.move_std(smoothed_relative, window=10, min_count=1).ravel() + smoothed_relative.ravel() - self.trpl_model([t], [g]).ravel()
        #fun_below = lambda g: 2*bn.move_std(smoothed_relative, window=10, min_count=1).ravel() - smoothed_relative.ravel() + self.trpl_model([t], [g]).ravel()
        #fun_above = lambda g: (2*std + smoothed_relative).ravel() - self.trpl_model([t], [g]).ravel()
        #fun_below = lambda g: (2*std - smoothed_relative).ravel() + self.trpl_model([t], [g]).ravel()

        #constraints.append({'type': 'ineq', 'fun': fun_above})
        constraints.append({'type': 'ineq', 'fun': fun_below})

        gamma_guess = self.fit_scaled_gamma([t], [smoothed_relative])[0][0]

        #import matplotlib.pyplot as plt
        #plt.plot((2*std + smoothed_relative).ravel())
        #plt.plot((-2*std + smoothed_relative).ravel())
        #plt.plot(self.trpl_model([t], [gamma_guess]).ravel())
        #plt.plot(1.1*smoothed_relative.ravel())
        #plt.plot(0.9*smoothed_relative.ravel())
        #plt.show()

        res = minimize(lambda g: -g, gamma_guess, constraints=constraints, options={'maxiter': 500})
        assert res.success

        print(gamma_guess)
        print(res.x[0])

        return res.x[0]
