import os
import re
from fpdf import FPDF

title = 'Malaria Model Validation Results'
section_title_size = 18
subsection_title_size = 11
body_text_size = 10
paragraph_spacing = 6
figure_width = 160


class PDF(FPDF):
    def header(self):
        # # Arial bold 15
        self.set_font('Arial', 'B', 15)
        # # Move to the right
        # self.cell(80)
        # # Title
        # # Arial bold 15
        # self.set_font('Arial', 'B', 25)
        # # Calculate width of title and position
        # w = self.get_string_width(title) + 6
        # self.set_x((210 - w) / 2)
        # # Colors of frame, background and text
        # self.set_draw_color(255, 255, 255)
        # self.set_fill_color(255, 255, 255)
        # self.set_text_color(44, 147, 194)
        # # self.set_text_color(60, 30, 240)
        # # Thickness of frame (1 mm)
        # self.set_line_width(1)
        # # Title
        # self.cell(w, 9, title, 1, 1, 'C', 1)
        # # Line break
        # self.ln(10)

    # Page footer
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')
        # Logo
        self.image('IDMlogo_small.png', 10, 282, 25)


# Instantiation of inherited class
pdf = PDF()
pdf.alias_nb_pages()
pdf.add_page()
# calculate page width
effective_page_width = pdf.w - 2 * pdf.l_margin

# _______  title  _________
# Move to the right
pdf.cell(80)
# Title
# Arial bold 15
pdf.set_font('Arial', 'B', 25)
# Calculate width of title and position
w = pdf.get_string_width(title) + 6
pdf.set_x((210 - w) / 2)
# Colors of frame, background and text
pdf.set_draw_color(255, 255, 255)
pdf.set_fill_color(255, 255, 255)
pdf.set_text_color(44, 147, 194)
# self.set_text_color(60, 30, 240)
# Thickness of frame (1 mm)
pdf.set_line_width(1)
# Title
pdf.cell(w, 9, title, 1, 1, 'C', 1)
# Line break
pdf.ln(10)

#________  overview of validation process  ____________
# pdf.set_font('Times', '', section_title_size)
# pdf.set_text_color(44, 147, 194)
# pdf.cell(0, 10, 'Validation overview', 0, 1)
# pdf.set_font('Times', '', body_text_size)
# pdf.set_text_color(0, 0, 0)
# pdf.multi_cell(0, 10, 'Insert very brief description of the validation process and where to look for more info.', 0, 1)
# # pdf.image('malaria_network_schematic.png', w=figure_width)
# pdf.ln(10)

#__________  brief summary  ________________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Results summary', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'Include some summary metrics on how the model performed for each of the validation relationships.', 0,
               1)
pdf.multi_cell(0, paragraph_spacing, 'PLACEHOLDER: Insert a table here.', 0, 1)
pdf.ln(10)

#______  age-incidence  _________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Incidence by age', 0, 1)
pdf.set_font('Times', '', subsection_title_size)
# matched-site simulations
pdf.cell(0, 10, 'Matched-site simulations', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below compare the age-incidence relationships from reference datasets and matched simulations.',
               0, 1)
pdf.image('_plots/site_compare_incidence_age.png', w=figure_width)
pdf.ln(10)
pdf.multi_cell(0, paragraph_spacing,
               'Next we look at the correlation between simulation and reference points (left plot) and at the correlation between the slopes of the lines between data points (right plot).',
               0, 1)
pdf.image('_plots/scatter_regression_incidence_age.png', w=figure_width)
pdf.ln(10)
pdf.multi_cell(0, paragraph_spacing,
               'The table below shows some summary metrics describing the match between the reference and simulation datasets.',
               0, 1)
pdf.multi_cell(0, paragraph_spacing, '{PLACEHOLDER: Insert a table with quantitative comparison results here.}', 0, 1)
# df = pd.read_csv('_plots/comparison_metric_table_incidence_age.csv')
# data = df.to_dict('records')
# line_height = body_text_size * 2.5
# col_width = effective_page_width / len(df.columns)  # distribute content evenly
# for row in data:
#     for k, v in row.items():
#         pdf.multi_cell(col_width, line_height, str(v), border=1, align="L", max_line_height=body_text_size)   # ln=3,
#     pdf.ln(line_height)
# for datum in row:
#     pdf.multi_cell(col_width, line_height, datum, border=1) #, new_x="RIGHT", new_y="TOP", max_line_height=body_text_size)
# pdf.ln(line_height)

# site sweeps
pdf.set_text_color(44, 147, 194)
pdf.set_font('Times', '', 15)
pdf.cell(0, 10, 'Sweep-site simulations', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below show a wider range of reference datasets plotted against a sweep of simulation-site characteristics.',
               0, 1)
pdf.multi_cell(0, paragraph_spacing, '{PLACEHOLDER: Insert sweep plots here.}', 0, 1)
pdf.ln(10)

#__________  age-prevalence  __________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Prevalence by age', 0, 1)
pdf.set_font('Times', '', subsection_title_size)
# matched-site simulations
pdf.cell(0, 10, 'Matched-site simulations', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below compare the age-prevalence relationships from reference datasets and matched simulations.',
               0, 1)
pdf.image('_plots/site_compare_prevalence_age.png', w=figure_width)
pdf.ln(10)
pdf.multi_cell(0, paragraph_spacing,
               'Next we look at the correlation between simulation and reference points (left plot) and at the correlation between the slopes of the lines between data points (right plot).',
               0, 1)
pdf.image('_plots/scatter_regression_prevalence_age.png', w=figure_width)
pdf.ln(10)
pdf.multi_cell(0, paragraph_spacing,
               'The table below shows some summary metrics describing the match between the reference and simulation datasets.',
               0, 1)
pdf.multi_cell(0, paragraph_spacing, '{PLACEHOLDER: Insert a table with quantitative comparison results here.}', 0, 1)

# site sweeps
pdf.set_text_color(44, 147, 194)
pdf.set_font('Times', '', 15)
pdf.cell(0, 10, 'Sweep-site simulations', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below show a wider range of reference datasets plotted against a sweep of simulation-site characteristics.',
               0, 1)
pdf.multi_cell(0, paragraph_spacing, '{PLACEHOLDER: Insert sweep plots here.}', 0, 1)
pdf.ln(10)

#____________  infectiousness  _________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Infectiousness to vectors', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, 10, 'Each plot panel corresponds to a site.', 0, 1)
files = [f for f in os.listdir('./_plots') if re.match(r'site_compare_infectiousness.*\.', f)]
for ff in files:
    site_name = ff.replace('.png', '')
    site_name = site_name.replace('site_compare_infectiousness_', '')
    pdf.multi_cell(0, 10, site_name, 0, 1)
    pdf.image('_plots/%s' % ff, w=figure_width * 3 / 4)
pdf.multi_cell(0, 10, 'Insert a table with quantitative comparison results here.', 0, 1)
pdf.ln(10)

#__________  infection duration  __________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Duration of infection', 0, 1)
pdf.set_font('Times', '', subsection_title_size)
# all ages
pdf.cell(0, 10, 'Across individuals of all ages', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below compare the duration over which individuals had positive tests in the reference dataset and matched simulations. The sampling design from the reference data was matched in the simulations.'
               'Observed infections are divided into two groups. "Censored" infections refer to infections where the individual was positive at the first or final survey of the study (so the infection may have extended beyond the period observed). '
               '"Start & finish observed" infections refer to infections were the individual was observed to have a negative test at the start and end of the infection.',
               0, 1)
pdf.image('_plots/site_compare_infect_duration_navrongo_2000.png', w=figure_width)
pdf.ln(10)
# all ages
pdf.cell(0, 10, 'Infection duration distribution by age', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, paragraph_spacing,
               'The plots below compare the duration over which individuals had positive tests in the reference dataset and matched simulations. The sampling design from the reference data was matched in the simulations.'
               'Observed infections are divided into two groups. "Censored" infections refer to infections where the individual was positive at the first or final survey of the study (so the infection may have extended beyond the period observed). '
               '"Start & finish observed" infections refer to infections were the individual was observed to have a negative test at the start and end of the infection.',
               0, 1)
pdf.image('_plots/site_compare_age_infect_duration_navrongo_2000.png', w=figure_width)
pdf.ln(10)

#__________  age-parasite density  ___________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Asexual parasite density by age', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, 10, 'Each plot panel corresponds to a site.', 0, 1)
files = [f for f in os.listdir('./_plots') if re.match(r'site_compare_par_dens_age.*\.', f)]
for ff in files:
    site_name = ff.replace('.png', '')
    site_name = site_name.replace('site_compare_par_dens_age_', '')
    pdf.multi_cell(0, 10, site_name, 0, 1)
    pdf.image('_plots/%s' % ff, w=figure_width * 1 / 2)
pdf.multi_cell(0, 10, 'Insert a table with quantitative comparison results here.', 0, 1)
pdf.ln(10)

#__________  age-gametocyte density  ___________
pdf.set_font('Times', '', section_title_size)
pdf.set_text_color(44, 147, 194)
pdf.cell(0, 10, 'Gametocyte density by age', 0, 1)
pdf.set_font('Times', '', body_text_size)
pdf.set_text_color(0, 0, 0)
pdf.multi_cell(0, 10, 'Each plot panel corresponds to a site.', 0, 1)
files = [f for f in os.listdir('./_plots') if re.match(r'site_compare_gamet_dens_age.*\.', f)]
for ff in files:
    site_name = ff.replace('.png', '')
    site_name = site_name.replace('site_compare_gamet_dens_age_', '')
    pdf.multi_cell(0, 10, site_name, 0, 1)
    pdf.image('_plots/%s' % ff, w=figure_width * 1 / 2)
pdf.multi_cell(0, 10, 'Insert a table with quantitative comparison results here.', 0, 1)
pdf.ln(10)

pdf.output('Malaria_model_validation_output2.pdf', 'F')
