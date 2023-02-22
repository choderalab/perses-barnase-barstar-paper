from pymbar import MBAR
from openmmtools.multistate import MultiStateSamplerAnalyzer, MultiStateReporter

class DataAnalyzer(object):

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            path to netcdf file containing simulation trajectory

        """
        self.reporter = MultiStateReporter(filename, open_mode='r')
        self.max_n_iterations = self.reporter.read_last_iteration()

    def get_equilibration_data(self, max_n_iterations=None):
        """
        Get equilibration data 

        Parameters
        ----------
        max_n_iterations : int, default None
            Maximum number of iterations to use when retrieving the equilibration data

        Returns
        -------
        equilibration_data : list
            contains number_equilibrated (number of iterations to discard due to equilibration),
            g_t (statistical inefficiency or subsample rate), 
            Neff_max (number of effective samples)
        
        """

        # Handle max_n_iterations argument
        if not max_n_iterations:
            max_n_iterations = self.max_n_iterations
        print(f"max_n_iterations: {max_n_iterations}")

        # Create analyzer
        analyzer = MultiStateSamplerAnalyzer(self.reporter, max_n_iterations=max_n_iterations)

        # Retrieve equilibration data
        equilibration_data = analyzer._get_equilibration_data()

        # Delete analyzer
        del analyzer

        return equilibration_data

    def get_free_energy(self, bootstrap=True, analyzer=None, max_n_iterations=None, n_equilibration_iterations=None, subsample_rate=None, initialize="zeros"):
        """
        Get free energy

        Parameters
        ----------
        bootstrap : boolean, default True
            whether to bootstrap the uncertainty
        analyzer : openmmtools.multistate.MultiStateSamplerAnalyzer, default None
            existing analyzer object to use for computing the free energy
            if None, will initialize a new one
        max_n_iterations : int, default None
            Maximum number of iterations to use when retrieving the equilibration data
        n_equilibration_iterations : int, default None
            number of iterations to discard due to equilibration
        subsample_rate : float, default None
            number of correlated samples in between decorrelated samples
        initialize : str, default "zeros"
            how to intialize the MBAR free energy estimate
            either with "zeros" or "BAR", the latter computes a BAR estimate to use for intializatio 

        Returns
        -------
        equilibration_data : list
            contains number_equilibrated (number of iterations until equilibration),
            g_t (statistical inefficiency or subsample rate),
            Neff_max (number of effective samples)

        """
        # If the analyzer was not supplied by the user
        if not analyzer:
             # Handle max_n_iterations argument
            if not max_n_iterations:
                max_n_iterations = self.max_n_iterations
            print(f"max_n_iterations: {max_n_iterations}")

            # Handle n_equilibration_iterations and subsample_rate arguments
            if (n_equilibration_iterations and not subsample_rate) or (not n_equilibration_iterations and subsample_rate):
                raise Exception("n_equilibration_iterations and subsample_rate must be specified together!")
            elif n_equilibration_iterations and subsample_rate:
                pass
            else:
                n_equilibration_iterations, subsample_rate, _ = self.get_equilibration_data(max_n_iterations=max_n_iterations)
            print(f"n_equilibration_iterations: {n_equilibration_iterations}")
            print(f"subsample rate: {subsample_rate}")

            # Create analyzer
            analyzer = MultiStateSamplerAnalyzer(self.reporter, max_n_iterations=max_n_iterations, n_equilibration_iterations=n_equilibration_iterations, statistical_inefficiency=subsample_rate)

        else:
            print("Analyzer supplied, ignoring user specified max_n_iterations, n_equilibration_iterations, and subsample_rate")

        if bootstrap:
            u_kn, N_k = analyzer._compute_mbar_decorrelated_energies()
            print(u_kn.shape)
            mbar = MBAR(u_kn, N_k, nbootstraps=200, solver_tolerance=1.0e-12, bootstrap_solver_tolerance=1.0e-6, initialize=initialize)
            Deltaf_ij, dDeltaf_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True, uncertainty_method='bootstrap')

        else: # If not bootstrapping, use solver_tolerance 1e-12
            u_kn, N_k = analyzer._compute_mbar_decorrelated_energies()
            print(u_kn.shape)
            mbar = MBAR(u_kn, N_k, solver_tolerance=1.0e-12, initialize=initialize)
            Deltaf_ij, dDeltaf_ij = mbar.getFreeEnergyDifferences(compute_uncertainty=True)

        # Create dict of results
        results = {"Deltaf": Deltaf_ij, "dDeltaf": dDeltaf_ij}

        # Create dict of metadata
        metadata = {"n_equilibration_iterations": n_equilibration_iterations,
                    "max_n_iterations": max_n_iterations,
                    "subsample_rate": subsample_rate,
                    "initialize": initialize}

        # Delete analyzer
        del analyzer

        return results, metadata
