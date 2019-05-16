__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"


# matplotlib.use('Agg')  # Not to use X server. For TravisCI.
import matplotlib.pyplot as plt  # noqa
import numpy as np
from pylatex import Document, Figure, NoEscape


def write_report(root, hdp, case_label):
    dof_label = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']

    width = r'1\textwidth'
    geometry_options = {"right": "2cm", "left": "2cm"}
    doc = Document(root.joinpath('report_{}'.format(case_label)), geometry_options=geometry_options)
    # print( hdp.available_dofs)
    # print(hdp.available_dirs)
    for sel_dof in hdp.available_dofs:
        with doc.create(Figure(position='htbp')) as plot:
            for sel_dir in hdp.available_dirs:
                plt.plot((2 * np.pi) / hdp.w, hdp.getRAO(sel_dof, sel_dir))
            plt.ylabel(dof_label[sel_dof - 1])

            plot.add_plot(width=NoEscape(width))
            plot.add_caption('RAO {}'.format(dof_label[sel_dof - 1]))
            plt.close()
    doc.generate_pdf(clean_tex=False)
