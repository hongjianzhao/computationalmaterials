from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected
from pymatgen.electronic_structure.core import Spin, OrbitalType

v = BSVasprun("vasprun.xml", parse_projected_eigen=True)
bs = v.get_band_structure(kpoints_filename="KPOINTS",line_mode=True)

emin=-8
emax=10

plt = BSPlotterProjected(bs)
plt.get_projected_plots_dots_patom_pmorb({'Zr':['dxy','dx2'],'O':['px','py','pz']}, {'Zr':[2],'O':[3]}, sum_atoms=None, sum_morbs=None, zero_to_efermi=True, ylim=[emin,emax], vbm_cbm_marker=False, selected_branches=None, w_h_size=(12, 8), num_column=None).show()
