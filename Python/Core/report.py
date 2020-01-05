__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"


# matplotlib.use('Agg')  # Not to use X server. For TravisCI.
import matplotlib.pyplot as plt  # noqa
from datetime import datetime
import numpy as np
from pylatex import Document, PageStyle, Head, Foot, MiniPage, \
    StandAloneGraphic, MultiColumn, Tabu, LongTabu, LargeText, MediumText, \
    LineBreak, NewPage, Tabularx, TextColor, simple_page_number, Command, \
    Figure, NoEscape, Section, LongTable
from pylatex.utils import bold, NoEscape

from pathlib import Path
import os

class DesignReport(Document):

    def __init__(self,settings):
        super().__init__()

        geometry_options = {
                "head": "40pt",

                "margin": "1.0in",
                "bottom": "1.0in",

                "includeheadfoot": True
        }
        self._doc = Document(settings._fio._case_dir.joinpath('report_{}'.format(settings._molo_model)),
                             geometry_options=geometry_options)

        # Generating first page style
        header = PageStyle("header")

        # Header image

        with header.create(Head("L")) as header_left:
            with header_left.create(MiniPage(width=NoEscape(r"0.49\textwidth"),
                                             pos='c')) as logo_wrapper:
                logo_file = str(settings.templates_dir.joinpath('MOLO_LOGO.png'))
                # logo_file= '{' + logo_file + '}'
                logo_file = logo_file.replace('\\', '/')
                # print(str(logo_file))
                logo_wrapper.append(StandAloneGraphic(image_options="width=120px",
                                                      filename=logo_file))
                # filename='logo.png'))

        # Add document title
        with header.create(Head("R")) as right_header:
            with right_header.create(MiniPage(width=NoEscape(r"0.49\textwidth"),
                                              pos='c', align='r')) as title_wrapper:
                title_wrapper.append('{}'.format(settings._molo_model))
                title_wrapper.append(LineBreak())
                title_wrapper.append(NoEscape(r'\today'))

        with header.create(Foot("L")):
            branch_address = MiniPage(
                    width=NoEscape(r"0.25\textwidth"),
                    pos='t', align='')
            branch_address.append(NoEscape(r"Verbun AS"))
            branch_address.append('\n')
            branch_address.append(NoEscape(r"R{\aa}dyrveien 32"))
            branch_address.append('\n')
            branch_address.append("1555 SON")
            branch_address.append('\n')
            branch_address.append("Norway")
            header.append(branch_address)

        with header.create(Foot("C")):
            today = datetime.today()
            header.append(NoEscape(r"Copyright {\textcopyright} 2019 Verbun AS"))


        with header.create(Foot("R")):
            header.append(simple_page_number())




        self._doc.preamble.append(header)
        self._doc.change_document_style("header")
        self._doc.add_color(name="lightgray", model="gray", description="0.80")
        self._doc.append(Command('author', 'Model: {}'.format(settings._molo_model)))
        self._doc.append(Command('title', 'MOLO Conceptual Core Report'))
        self._doc.append(Command('maketitle'))
        self._doc.append(NewPage())
        self._doc.append(Command('tableofcontents'))

        # Add Heading





        # with self._doc.create(MiniPage(align='c')):
        #     self._doc.append(LargeText(bold("MOLO Conceptual Core Report")))
        #     self._doc.append(LineBreak())
        #     self._doc.append(MediumText(bold('{}'.format(settings._molo_label))))



    def write_hydrostatic_report_latex_table(self, hydrostatics):
        def build_table_data(text, data, precision=3, dtype='f'):
            textwidth = 40
            label_template = "{:{fill}{align}{width}}"
            num_template = "{:{align}{width}.{precision}{dtype}} "
            strings = list()
            # strings.append(label_template.format(text, fill=' ', align='<', width=40))
            sa = label_template.format(text, fill=' ', align='<', width=40)
            sb = 'No data'
            # for i in range(data.shape()[0]):
            if type(data) is np.ndarray:
                data = data.tolist()
                for val in data:
                    sb = num_template.format(val, align='<', width=8, precision=precision, dtype=dtype)

            else:
                sb = num_template.format(data, align='<', width=8, precision=precision, dtype=dtype)
            return [sa, sb]
        doc=self._doc
        doc.append(NewPage())
        with doc.create(Section('Hydrostatic report')):
            # Generate hydrostatic data table
            with doc.create(LongTable("l l")) as data_table:
                data_table.add_row(build_table_data('Gravity acceleration (M/S**2)', hydrostatics.gravity, precision=2))
                data_table.add_row(build_table_data('Density of water (kg/M**3)', hydrostatics.rho_water, precision=1))

                data_table.add_row(build_table_data('Waterplane area (M**2)', hydrostatics.flotation_surface_area, precision=1))
                data_table.add_row(build_table_data('Waterplane center (M)', hydrostatics.flotation_center[:2], precision=3))
                data_table.add_row(build_table_data('Wet area (M**2)', hydrostatics.wet_surface_area, precision=1))
                data_table.add_row(build_table_data('Displacement volume (M**3)', hydrostatics.displacement_volume, precision=3))
                data_table.add_row(build_table_data('Displacement mass (tons)', hydrostatics.displacement, precision=3))
                data_table.add_row(build_table_data('Buoyancy center (M)', hydrostatics.buoyancy_center, precision=3))
                data_table.add_row(build_table_data('Center of gravity (M)', hydrostatics.gravity_center, precision=3))

                data_table.add_row(build_table_data('Draught (M)', hydrostatics.hs_data['draught'], precision=3))
                data_table.add_row(build_table_data('Length overall submerged (M)', hydrostatics.hs_data['los'], precision=2))
                data_table.add_row(build_table_data('Breadth overall submerged (M)', hydrostatics.hs_data['bos'], precision=2))
                data_table.add_row(build_table_data('Length at Waterline LWL (M)', hydrostatics.hs_data['lwl'], precision=2))
                data_table.add_row(build_table_data('Forward perpendicular (M)', hydrostatics.hs_data['fp'], precision=2))

                data_table.add_row(['',''])
                data_table.add_row(
                    build_table_data('Transversal metacentric radius (M)', hydrostatics.transversal_metacentric_radius, precision=3))
                data_table.add_row(
                    build_table_data('Transversal metacentric height GMt (M)', hydrostatics.transversal_metacentric_height,
                               precision=3))
                data_table.add_row(
                    build_table_data('Longitudinal metacentric radius (M)', hydrostatics.longitudinal_metacentric_radius, precision=3))
                data_table.add_row(
                    build_table_data('Longitudinal metacentric height GMl (M)', hydrostatics.longitudinal_metacentric_height,
                               precision=3))


                data_table.add_row(['',''])
                data_table.add_row(['HYDROSTATIC STIFFNESS COEFFICIENTS:',''])
                data_table.add_row(build_table_data('K33 (N/M)', hydrostatics.S33, precision=4, dtype='E'))
                data_table.add_row(build_table_data('K34 (N)', hydrostatics.S34, precision=4, dtype='E'))
                data_table.add_row(build_table_data('K35 (N)', hydrostatics.S35, precision=4, dtype='E'))
                data_table.add_row(build_table_data('K44 (N.M)', hydrostatics.S44, precision=4, dtype='E'))
                data_table.add_row(build_table_data('K45 (N.M)', hydrostatics.S45, precision=4, dtype='E'))
                data_table.add_row(build_table_data('K55 (N.M)', hydrostatics.S55, precision=4, dtype='E'))

                data_table.add_row(['', ''])
                data_table.add_row(['INERTIAS:', '']) #TODO: Use intertias from mass model instead
                #data_table.add_row('\tINERTIAS:\n')
                data_table.add_row(build_table_data('Ixx', hydrostatics.hs_data['Ixx'], precision=3, dtype='E'))
                data_table.add_row(build_table_data('Ixy', hydrostatics.hs_data['Ixy'], precision=3, dtype='E'))
                data_table.add_row(build_table_data('Ixz', hydrostatics.hs_data['Ixz'], precision=3, dtype='E'))
                data_table.add_row(build_table_data('Iyy', hydrostatics.hs_data['Iyy'], precision=3, dtype='E'))
                data_table.add_row(build_table_data('Iyz', hydrostatics.hs_data['Iyz'], precision=3, dtype='E'))
                data_table.add_row(build_table_data('Izz', hydrostatics.hs_data['Izz'], precision=3, dtype='E'))

                data_table.add_row(['', ''])
                data_table.add_row(['RESIDUALS:', ''])
                data_table.add_row(build_table_data('Absolute', hydrostatics.residual, precision=3, dtype='E'))
                data_table.add_row(build_table_data('Relative', hydrostatics.residual / hydrostatics._scale, precision=3, dtype='E'))
                data_table.add_row(build_table_data('Relative tolerance', hydrostatics.reltol, precision=1, dtype='E'))


    def write_hydrostatic_report_latex_list(self, hydrostatics):
        """Returns a hydrostatic report for the current configuration

         Returns
         -------
         str
         """

        def hspace():
            return '\n'

        # def build_line(text, data, precision=3, dtype='f'):
        #     textwidth = 40
        #     try:
        #         line = '\t{:-<{textwidth}}>  {:< .{precision}{dtype}}\n'.format(str(text).upper(), data,
        #                                                                         precision=precision,
        #                                                                         textwidth=textwidth,
        #                                                                         dtype=dtype
        #                                                                         )
        #     except ValueError:
        #         if isinstance(data, np.ndarray):
        #             if data.ndim == 1:
        #                 data_str = ''.join(['{:< 10.{precision}{dtype}}'.format(val, precision=precision, dtype=dtype)
        #                                     for val in data])
        #                 line = '\t{:-<{textwidth}}>  {}\n'.format(str(text).upper(), data_str, textwidth=textwidth)
        #
        #                 # else:
        #                 #     print data
        #
        #     return line

        def build_line(text, data, precision=3, dtype='f'):
            textwidth = 40
            label_template = "{:{fill}{align}{width}}"
            num_template = "{:{align}{width}.{precision}{dtype}} "
            strings = list()
            strings.append(label_template.format(text, fill=' ', align='<', width=40))
            # for i in range(data.shape()[0]):
            if type(data) is np.ndarray:
                data = data.tolist()
                for val in data:
                    strings.append(num_template.format(val, align='<', width=8, precision=precision, dtype=dtype))

            else:
                strings.append(num_template.format(data, align='<', width=8, precision=precision, dtype=dtype))
            return ''.join(strings) + '\n'


        msg = '\n'
        # title = 'Hydrostatic report ({0})\n
        # \tGenerated by meshmagick on {1}'.format(self.mesh.name, strftime('%c')).upper()
        # msg += header(title)

        msg += hspace()
        msg += build_line('Gravity acceleration (M/S**2)', self.gravity, precision=2)
        msg += build_line('Density of water (kg/M**3)', self.rho_water, precision=1)

        msg += hspace()
        msg += build_line('Waterplane area (M**2)', self.flotation_surface_area, precision=1)
        msg += build_line('Waterplane center (M)', self.flotation_center[:2], precision=3)
        msg += build_line('Wet area (M**2)', self.wet_surface_area, precision=1)
        msg += build_line('Displacement volume (M**3)', self.displacement_volume, precision=3)
        msg += build_line('Displacement mass (tons)', self.displacement, precision=3)
        msg += build_line('Buoyancy center (M)', self.buoyancy_center, precision=3)
        msg += build_line('Center of gravity (M)', self.gravity_center, precision=3)

        msg += hspace()
        msg += build_line('Draught (M)', self.hs_data['draught'], precision=3)
        msg += build_line('Length overall submerged (M)', self.hs_data['los'], precision=2)
        msg += build_line('Breadth overall submerged (M)', self.hs_data['bos'], precision=2)
        msg += build_line('Length at Waterline LWL (M)', self.hs_data['lwl'], precision=2)
        msg += build_line('Forward perpendicular (M)', self.hs_data['fp'], precision=2)

        msg += hspace()
        msg += build_line('Transversal metacentric radius (M)', self.transversal_metacentric_radius, precision=3)
        msg += build_line('Transversal metacentric height GMt (M)', self.transversal_metacentric_height, precision=3)
        msg += build_line('Longitudinal metacentric radius (M)', self.longitudinal_metacentric_radius, precision=3)
        msg += build_line('Longitudinal metacentric height GMl (M)', self.longitudinal_metacentric_height, precision=3)

        msg += hspace()
        msg += '\tHYDROSTATIC STIFFNESS COEFFICIENTS:\n'
        msg += build_line('K33 (N/M)', self.S33, precision=4, dtype='E')
        msg += build_line('K34 (N)', self.S34, precision=4, dtype='E')
        msg += build_line('K35 (N)', self.S35, precision=4, dtype='E')
        msg += build_line('K44 (N.M)', self.S44, precision=4, dtype='E')
        msg += build_line('K45 (N.M)', self.S45, precision=4, dtype='E')
        msg += build_line('K55 (N.M)', self.S55, precision=4, dtype='E')

        # Il faut faire une correction avec le plan de la flottaison de certains coeffs
        msg += hspace()
        msg += '\tINERTIAS:\n'
        msg += build_line('Ixx', self.hs_data['Ixx'], precision=3, dtype='E')
        msg += build_line('Ixy', self.hs_data['Ixy'], precision=3, dtype='E')
        msg += build_line('Ixz', self.hs_data['Ixz'], precision=3, dtype='E')
        msg += build_line('Iyy', self.hs_data['Iyy'], precision=3, dtype='E')
        msg += build_line('Iyz', self.hs_data['Iyz'], precision=3, dtype='E')
        msg += build_line('Izz', self.hs_data['Izz'], precision=3, dtype='E')

        msg += hspace()
        msg += '\tRESIDUALS:\n'
        msg += build_line('Absolute', self.residual, precision=3, dtype='E')
        msg += build_line('Relative', self.residual / self._scale, precision=3, dtype='E')
        msg += build_line('Relative tolerance', self.reltol, precision=1, dtype='E')
        # residual = self.residual
        # msg += ('\nResidual:\n')
        # msg += ('Delta Fz = %.3f N\n' % residual[0])
        # msg += ('Delta Mx = %.3f Nm\n' % residual[1])
        # msg += ('Delta My = %.3f Nm\n' % residual[2])
        #
        # rel_res = residual / self._scale
        # msg += ('\nRelative residual:\n')
        # msg += ('Delta Fz = %E\n' % rel_res[0])
        # msg += ('Delta Mx = %E\n' % rel_res[1])
        # msg += ('Delta My = %E\n' % rel_res[2])
        # msg += ('Relative tolerance of the solver: %.1E\n' % self.reltol)

        return msg


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
                plt.plot((2 * np.pi) / hdp.w, hdp.calc_rao_linear(sel_dof, sel_dir))
            plt.ylabel(dof_label[sel_dof - 1])

            plot.add_plot(width=NoEscape(width))
            plot.add_caption('RAO {}'.format(dof_label[sel_dof - 1]))
            plt.close()
    doc.generate_pdf(clean_tex=False)



