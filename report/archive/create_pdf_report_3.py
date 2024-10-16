import os
import re
import pandas as pd
from fpdf import FPDF, TitleStyle, XPos, YPos
from datetime import datetime
import simulations.manifest as manifest

now = datetime.now()
now_str = now.strftime("%c")
date, time = now.strftime("%m-%d-%Y %H-%M-%S").split(' ')

title = 'Malaria Model Validation Results'
pdf_name = f'Malaria_model_validation_output_{date}_{time}.pdf'
plot_folder = '_plots_R'


def new_section(pdf: FPDF, text: str, image_list: list = None, table_name: str = None, new_x: int = XPos.LMARGIN,
                new_y: int = YPos.NEXT, **kwargs):
    pdf.multi_cell(
        w=pdf.epw,
        h=pdf.font_size,
        txt=text,
        new_x=new_x,
        new_y=new_y,
        **kwargs,
    )
    if image_list:
        pdf.ln(10)
        for image_name in image_list:
            pdf.image(image_name, w=185)
    if table_name:
        pdf.ln(10)
        df = pd.read_csv(table_name)
        line_height = 8
        col_width = pdf.epw / len(df.columns)  # distribute content evenly
        pdf.set_font(size=8)
        # header
        pdf.set_font(style="B")  # enabling bold text
        for col_name in df.columns:
            pdf.cell(col_width, line_height, col_name, border=1)
        pdf.ln(line_height)
        pdf.set_font(style="")  # disabling bold text
        for index, row in df.iterrows():
            for datum in row:
                if isinstance(datum, (int, float)):
                    datum = str(round(datum, 2))
                pdf.multi_cell(col_width, line_height, datum, border=1,
                               new_x=XPos.RIGHT, new_y=YPos.TOP, max_line_height=pdf.font_size)
            pdf.ln(line_height)
        pdf.set_font(size=12)


def render_toc(pdf, outline):
    pdf.x = 10
    pdf.y += 20
    pdf.set_font("Helvetica", size=16)
    pdf.underline = True
    new_section(pdf, "Table of contents")
    pdf.underline = False
    pdf.y += 5
    pdf.set_font("Courier", size=10)
    for section in outline:
        link = pdf.add_link()
        pdf.set_link(link, page=section.page_number)
        text = f'{" " * section.level * 2} {section.name}'
        text += (
            f' {"." * (70 - section.level * 2 - len(section.name))} {section.page_number}'
        )
        pdf.multi_cell(
            w=pdf.epw,
            h=pdf.font_size,
            txt=text,
            new_x=XPos.LMARGIN,
            new_y=YPos.NEXT,
            align="L",
            link=link,
        )


class PDF(FPDF):
    def header(self):
        # Logo
        self.image('IDMlogo_small.png', 10, 8, 33)
        # Arial bold 15
        self.set_font('Arial', 'B', 15)
        # Move to the right
        self.cell(80)
        # Title
        # Arial bold 15
        self.set_font('Arial', 'B', 12)
        # Calculate width of title and position
        w = self.get_string_width(title) + 6
        self.set_x((210 - w) / 2)
        # Colors of frame, background and text
        self.set_draw_color(255, 255, 255)
        self.set_fill_color(255, 255, 255)
        self.set_text_color(44, 147, 194)
        # Thickness of frame (1 mm)
        self.set_line_width(1)
        # Title
        self.cell(w, txt=title)
        # Line break
        self.ln(20)

    # Page footer
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, txt='Page ' + str(self.page_no()) + '/{nb}' + f'\t\t{now_str}')


# Instantiation of inherited class
pdf = PDF()
pdf.alias_nb_pages()

# add outlines
pdf.set_font("Helvetica")
pdf.set_section_title_styles(
    # Level 0 titles:
    TitleStyle(
        font_family="Helvetica",
        font_style="B",
        font_size_pt=16,
        color=(44, 147, 194),
        underline=True,
        t_margin=10,
        l_margin=10,
        b_margin=0,
    ),
    # Level 1 subtitles:
    TitleStyle(
        font_family="Helvetica",
        font_style="B",
        font_size_pt=14,
        color=(44, 147, 194),
        underline=True,
        t_margin=10,
        l_margin=10,
        b_margin=5,
    ),
)
pdf.add_page()
pdf.set_y(50)
pdf.set_font(size=25)
pdf.set_text_color(44, 147, 194)
new_section(pdf, 'Validation Report', align="C")
pdf.set_font(size=16)

suite_id_filepath = manifest.CURRENT_DIR / manifest.suite_id_file
if suite_id_filepath.is_file():
    with open(suite_id_filepath, 'r') as suite_id_file:
        suite_id = suite_id_file.readline()
else:
    suite_id = 'NA'

version_file_filepath = manifest.CURRENT_DIR / manifest.version_file
if version_file_filepath.is_file():
    with open(version_file_filepath, 'r') as version_file:
        era_version = version_file.readline().rstrip()
        era_branch = version_file.readline().rstrip()
        emodpy_malaria_version = version_file.readline().rstrip()

else:
    era_version = 'Eradication version: ___'
    era_branch = 'branch: ___'
    emodpy_malaria_version = 'emodpy-malaria version: ___'

new_section(pdf, f'\n {date} \n '
                 f'{era_version} \n '
                 f'{era_branch} \n '
                 f'{emodpy_malaria_version} \n '
                 f'Suite ID: {suite_id}', align="C")
pdf.set_font(size=12)
pdf.set_text_color(0, 0, 0)
pdf.insert_toc_placeholder(render_toc)

# define sections and content
section_and_content = {'Introduction':
    {'Background':
        [
            'The goal of this report is to help users quickly identify whether updated versions of '
            'the malaria model are still well-calibrated to capture a range of relevant real-world '
            'malaria observations. \n\n'
            'The figures and tables compare simulation output generated with a particular version of '
            'the Eradication.exe and of emodpy-malaria with 1) the simulation results generated by '
            'earlier versions of Eradication.exe and emodpy-malaria (the versions used to calibrate '
            'the model) and 2) reference datasets from real-world observations. \n\n'
            'This report was generated by running the malaria model validation workflow available at '
            'https://github.com/InstituteforDiseaseModeling/malaria-model_validation. '
            'Additional information on the reference datasets and on the simulation assumptions '
            'are available from the repo in "Notes on reference datasets and simulation '
            'assumptions.docx," and instructions on how to re-run the validation comparisons are '
            'in the README file.',
            None,
            None]},
    'Results summary':
        {'Performance compared to model version from calibration':
            [
                'The table below shows, for each validation relationship (rows), the mean absolute difference between '
                'all reference and matched simulation datapoints for both the new (first column) and benchmark (second '
                'column) simulations. The final three columns of the table show the number of sites where the new '
                'simulations matched the reference dataset better, similarly, or worse compared to the benchmark '
                'simulations.'
                ''
                '\n\n',
                None,
                f'{plot_folder}/summary_table_sim_benchmark.csv'],
         },
    'Visual comparison of reference data and matched simulations':
        {'Incidence by age':
            [
                'The plots below compare the age-incidence relationships from reference datasets and matched simulations.',
                [f'{plot_folder}/site_compare_incidence_age.png'],
                None],
            'Prevalence by age':
                [
                    'The plots below compare the age-prevalence relationships from reference datasets and matched simulations.',
                    [f'{plot_folder}/site_compare_prevalence_age.png'],
                    None],
            'Infectiousness to vectors':
                [
                    'Each of the below plot panels corresponds to a site. Within a plot panel, each row corresponds to an age group and each column corresponds to the month when sampling occurred.\n'
                    'The x-axis shows the gametocyte density in an infection. The y-axis shows how infectious an individual is to mosquitoes. '
                    'The dot size shows how often a person of a given age and gametocyte density falls into each of the infectiousness bins '
                    "(each column's dot sizes sum to one).\n"
                    "In the reference datasets, the sample size is sometimes quite small.",
                    [f'{plot_folder}/%s' % ff for ff in os.listdir(f'./{plot_folder}') if
                     re.match(r'site_compare_infectiousness.*\.', ff)],
                    None],
            'Duration of infection - all ages':
                [
                    'The plots below compare the duration over which individuals had positive tests in the reference dataset and matched simulations. The sampling design from the reference data was matched in the simulations.'
                    'Observed infections are divided into two groups. "Censored" infections refer to infections where the individual was positive at the first or final survey of the study (so the infection may have extended beyond the period observed). '
                    '"Start & finish observed" infections refer to infections were the individual was observed to have a negative test at the start and end of the infection. The two types of infection duration records are illustrated in the figure below.',
                    ['infection_duration_censoring_illustration.png',
                     f'{plot_folder}/site_compare_infect_duration_navrongo_2000.png'],
                    None],
            'Duration of infection - by age':
                [
                    'The plots below compare the duration over which individuals had positive tests in the reference dataset and matched simulations. The sampling design from the reference data was matched in the simulations.'
                    'Observed infections are divided into two groups. "Censored" infections refer to infections where the individual was positive at the first or final survey of the study (so the infection may have extended beyond the period observed). '
                    '"Start & finish observed" infections refer to infections were the individual was observed to have a negative test at the start and end of the infection. The two types of infection duration records are illustrated in the figure below.\n'
                    'In the plot panel below, columns correspond to the age group (in years) and rows correspond to whether or not the start and end of the infection was observed.',
                    [f'{plot_folder}/site_compare_infect_duration_age_navrongo_2000.png'],
                    None],
            'Asexual parasite density by age':
                ['The plots below compare the distribution of parasite densities across ages and '
                 'seasons from reference datasets and matched simulations. Each plot panel corresponds '
                 'to a site. Note that some of the reference datasets have small sample sizes, '
                 'especially in the youngest age groups.',
                 [f'{plot_folder}/%s' % ff for ff in os.listdir(f'./{plot_folder}') if
                  re.match(r'site_compare_asex_dens_age.*\.', ff)],
                 None],
            'Gametocyte density by age':
                ['The plots below compare the distribution of gametocyte densities across ages and '
                 'seasons from reference datasets and matched simulations. Each plot panel corresponds '
                 'to a site. Note that some of the reference datasets have small sample sizes, '
                 'especially in the youngest age groups.',
                 [f'{plot_folder}/%s' % ff for ff in os.listdir(f'./{plot_folder}') if
                  re.match(r'site_compare_gamet_dens_age_.*\.', ff)],
                 None]},
    'Additional comparisons of reference data and matched simulations':
        {'Incidence by age':
             ['Below, the plots show the correlation between the simulation and reference incidence '
              'in each age group (left plot) and the correlation between the simulation and '
              'reference slopes when moving between the incidence in one age group and the incidence '
              'in the next oldest age group (right plot).',
              [f'{plot_folder}/scatter_regression_incidence_age.png'],
              f'{plot_folder}/comparison_metric_table_incidence_age.csv'],
         'Prevalence by age':
             [
                 'Below, the plots show the correlation between the simulation and reference prevalence '
                 'in each age group (left plot) and the correlation between the simulation and '
                 'reference slopes when moving between the prevalence in one age group and the prevalence '
                 'in the next oldest age group (right plot).',
                 [f'{plot_folder}/scatter_regression_prevalence_age.png'],
                 f'{plot_folder}/comparison_metric_table_prevalence_age.csv'],
         # 'Infectiousness to vectors':
         #     ['...',
         #      None,
         #      None],
         # 'Duration of infection':
         #     ['...',
         #      None,
         #      None],
         # 'Asexual parasite density by age':
         #     ['...',
         #      None,
         #      None],
         # 'Gametocyte density by age':
         #     ['...',
         #      None,
         #      None]
         },
    # 'Visual comparison of reference data and simulation sweeps':
    #     {'':
    #          ['...\n...',
    #           None,
    #           None],
    #      '':
    #          ['...',
    #           None,
    #           None]},

    'Comparisons from prior EMOD publications':
        {'Incidence and prevalence by age':
             ['The top plots come from McCarthy et al. 2015 and the bottom plots come from Selvaraj et al. 2018.',
              ['_prior_recalibration_published_figures/McCarthy_2015_age_pfpr_incidence.png',
               '_prior_recalibration_published_figures/Selvaraj_2018_age_incidence.png'],
              None],
         'Infectiousness to vectors':
             ['The following plot comes from Selvaraj et al. 2018',
             ['_prior_recalibration_published_figures/Selvaraj_2018_infectiousness.png'],
             None],
         'Parasite densities':
             ['The following plots come from Selvaraj et al. 2018',
             ['_prior_recalibration_published_figures/Selvaraj_2018_parasite_densities.png',
             '_prior_recalibration_published_figures/Selvaraj_2018_parasite_densities2.png',
             '_prior_recalibration_published_figures/Selvaraj_2018_parasite_densities3.png'],
             None]},
}

for i, (section_tile, section_content) in enumerate(section_and_content.items()):
    pdf.start_section(f"{i + 1}. {section_tile}")
    subsection_names, text_and_images = section_content.keys(), section_content.values()
    for j, subsection_name, text_and_images in zip(range(len(subsection_names)), subsection_names, text_and_images):
        section_text, image_list, table_name = text_and_images
        pdf.start_section(f'{i + 1}.{j + 1} ' + subsection_name, level=1)
        new_section(
            pdf,
            section_text,
            image_list=image_list,
            table_name=table_name
        )
        pdf.add_page()
new_section(pdf, 'This page left unintentionally blank', align="C")

# pdf.add_page()
# for section_name, section_content, image_name in zip(section_names, section_content_list, image_names):
#     pdf.new_section(section_name=section_name,
#                     section_content=section_content,
#                     image_name=image_name)

pdf.output(pdf_name)
