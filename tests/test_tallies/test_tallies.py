#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)

from testing_harness import PyAPITestHarness
from openmc.filter import *
from openmc import Mesh, Tally, Tallies
from openmc.source import Source
from openmc.stats import Box


class TalliesTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Build default materials/geometry
        self._input_set.build_default_materials_and_geometry()

        # Set settings explicitly
        self._input_set.settings.batches = 5
        self._input_set.settings.inactive = 0
        self._input_set.settings.particles = 400
        self._input_set.settings.source = Source(space=Box(
            [-160, -160, -183], [160, 160, 183]))

        azimuthal_bins = (-3.14159, -1.8850, -0.6283, 0.6283, 1.8850, 3.14159)
        azimuthal_filter = AzimuthalFilter(azimuthal_bins)
        azimuthal_tally1 = Tally()
        azimuthal_tally1.filters = [azimuthal_filter]
        azimuthal_tally1.scores = ['flux']
        azimuthal_tally1.estimator = 'tracklength'

        azimuthal_tally2 = Tally()
        azimuthal_tally2.filters = [azimuthal_filter]
        azimuthal_tally2.scores = ['flux']
        azimuthal_tally2.estimator = 'analog'

        mesh_2x2 = Mesh(mesh_id=1)
        mesh_2x2.lower_left  = [-182.07, -182.07]
        mesh_2x2.upper_right = [182.07,  182.07]
        mesh_2x2.dimension = [2, 2]
        mesh_filter = MeshFilter(mesh_2x2)
        azimuthal_tally3 = Tally()
        azimuthal_tally3.filters = [azimuthal_filter, mesh_filter]
        azimuthal_tally3.scores = ['flux']
        azimuthal_tally3.estimator = 'tracklength'

        cellborn_tally = Tally()
        cellborn_tally.filters = [CellbornFilter((10, 21, 22, 23))]
        cellborn_tally.scores = ['total']

        dg_tally = Tally()
        dg_tally.filters = [DelayedGroupFilter((1, 2, 3, 4, 5, 6))]
        dg_tally.scores = ['delayed-nu-fission']

        four_groups = (0.0, 0.253, 1.0e3, 1.0e6, 20.0e6)
        energy_filter = EnergyFilter(four_groups)
        energy_tally = Tally()
        energy_tally.filters = [energy_filter]
        energy_tally.scores = ['total']

        energyout_filter = EnergyoutFilter(four_groups)
        energyout_tally = Tally()
        energyout_tally.filters = [energyout_filter]
        energyout_tally.scores = ['scatter']

        transfer_tally = Tally()
        transfer_tally.filters = [energy_filter, energyout_filter]
        transfer_tally.scores = ['scatter', 'nu-fission']

        material_tally = Tally()
        material_tally.filters = [MaterialFilter((1, 2, 3, 4))]
        material_tally.scores = ['total']

        mu_bins = (-1.0, -0.5, 0.0, 0.5, 1.0)
        mu_filter = MuFilter(mu_bins)
        mu_tally1 = Tally()
        mu_tally1.filters = [mu_filter]
        mu_tally1.scores = ['scatter', 'nu-scatter']

        mu_tally2 = Tally()
        mu_tally2.filters = [mu_filter, mesh_filter]
        mu_tally2.scores = ['scatter', 'nu-scatter']

        polar_bins = (0.0, 0.6283, 1.2566, 1.8850, 2.5132, 3.14159)
        polar_filter = PolarFilter(polar_bins)
        polar_tally1 = Tally()
        polar_tally1.filters = [polar_filter]
        polar_tally1.scores = ['flux']
        polar_tally1.estimator = 'tracklength'

        polar_tally2 = Tally()
        polar_tally2.filters = [polar_filter]
        polar_tally2.scores = ['flux']
        polar_tally2.estimator = 'analog'

        polar_tally3 = Tally()
        polar_tally3.filters = [polar_filter, mesh_filter]
        polar_tally3.scores = ['flux']
        polar_tally3.estimator = 'tracklength'

        universe_tally = Tally()
        universe_tally.filters = [UniverseFilter((1, 2, 3, 4, 6, 8))]
        universe_tally.scores = ['total']

        cell_filter = CellFilter((10, 21, 22, 23, 60))
        score_tallies = [Tally(), Tally(), Tally()]
        for t in score_tallies:
            t.filters = [cell_filter]
            t.scores = ['absorption', 'delayed-nu-fission', 'events', 'fission',
                        'inverse-velocity', 'kappa-fission', '(n,2n)', '(n,n1)',
                        '(n,gamma)', 'nu-fission', 'scatter', 'elastic',
                        'total', 'prompt-nu-fission', 'fission-q-prompt',
                        'fission-q-recoverable']
        score_tallies[0].estimator = 'tracklength'
        score_tallies[1].estimator = 'analog'
        score_tallies[2].estimator = 'collision'

        cell_filter2 = CellFilter((21, 22, 23, 27, 28, 29, 60))
        flux_tallies = [Tally() for i in range(4)]
        for t in flux_tallies:
            t.filters = [cell_filter2]
        flux_tallies[0].scores = ['flux']
        for t in flux_tallies[1:]:
            t.scores = ['flux-y5']
        flux_tallies[1].estimator = 'tracklength'
        flux_tallies[2].estimator = 'analog'
        flux_tallies[3].estimator = 'collision'

        scatter_tally1 = Tally()
        scatter_tally1.filters = [cell_filter]
        scatter_tally1.scores = ['scatter', 'scatter-1', 'scatter-2', 'scatter-3',
                                 'scatter-4', 'nu-scatter', 'nu-scatter-1',
                                 'nu-scatter-2', 'nu-scatter-3', 'nu-scatter-4']

        scatter_tally2 = Tally()
        scatter_tally2.filters = [cell_filter]
        scatter_tally2.scores = ['scatter-p4', 'scatter-y4', 'nu-scatter-p4',
                                 'nu-scatter-y3']

        total_tallies = [Tally() for i in range(4)]
        for t in total_tallies:
            t.filters = [cell_filter]
        total_tallies[0].scores = ['total']
        for t in total_tallies[1:]:
            t.scores = ['total-y4']
            t.nuclides = ['U235', 'total']
        total_tallies[1].estimator = 'tracklength'
        total_tallies[2].estimator = 'analog'
        total_tallies[3].estimator = 'collision'

        all_nuclide_tallies = [Tally(), Tally()]
        for t in all_nuclide_tallies:
            t.filters = [cell_filter]
            t.nuclides = ['all']
            t.scores = ['total']
        all_nuclide_tallies[0].estimator = 'tracklength'
        all_nuclide_tallies[0].estimator = 'collision'

        self._input_set.tallies = Tallies()
        self._input_set.tallies += (
            [azimuthal_tally1, azimuthal_tally2, azimuthal_tally3,
             cellborn_tally, dg_tally, energy_tally, energyout_tally,
             transfer_tally, material_tally, mu_tally1, mu_tally2,
             polar_tally1, polar_tally2, polar_tally3, universe_tally])
        self._input_set.tallies += score_tallies
        self._input_set.tallies += flux_tallies
        self._input_set.tallies += (scatter_tally1, scatter_tally2)
        self._input_set.tallies += total_tallies
        self._input_set.tallies += all_nuclide_tallies

        self._input_set.export()

    def _get_results(self):
        return super(TalliesTestHarness, self)._get_results(hash_output=True)

    def _cleanup(self):
        super(TalliesTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = TalliesTestHarness('statepoint.5.*', True)
    harness.main()
