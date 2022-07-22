from fpdf import FPDF, TitleStyle, XPos, YPos
from datetime import datetime

now = datetime.now()
now_str = now.strftime("%c")
date, time = now.strftime("%d-%m-%Y %H-%M-%S").split(' ')

title = 'Malaria Model Doc'
pdf_name = f'Malaria_model_validation_output_{date}_{time}.pdf'


def new_section(pdf: FPDF, text: str, image_name: str = None, new_x: int = XPos.LMARGIN, new_y: int = YPos.NEXT,
                **kwargs):
    pdf.multi_cell(
        w=pdf.epw,
        h=pdf.font_size,
        txt=text,
        new_x=new_x,
        new_y=new_y,
        **kwargs,
    )
    if image_name:
        pdf.image(image_name, w=150)


def render_toc(pdf, outline):
    pdf.x = 10
    pdf.y += 50
    pdf.set_font("Helvetica", size=16)
    pdf.underline = True
    new_section(pdf, "Table of contents:")
    pdf.underline = False
    pdf.y += 20
    pdf.set_font("Courier", size=12)
    for section in outline:
        link = pdf.add_link()
        pdf.set_link(link, page=section.page_number)
        text = f'{" " * section.level * 2} {section.name}'
        text += (
            f' {"." * (60 - section.level * 2 - len(section.name))} {section.page_number}'
        )
        pdf.multi_cell(
            w=pdf.epw,
            h=pdf.font_size,
            txt=text,
            new_x=XPos.LMARGIN,
            new_y=YPos.NEXT,
            align="C",
            link=link,
        )


class PDF(FPDF):
    def header(self):
        # Logo
        self.image('logo_IDM_horz3_RGB-min.jpg', 10, 8, 33)
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
        self.set_text_color(60, 30, 240)
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

    def new_section(self, section_name, section_content, image_name):
        pdf.set_font('Times', '', 24)
        pdf.cell(0, 30, txt=section_name, new_x=XPos.LEFT, new_y=YPos.NEXT)
        pdf.set_font('Times', '', 12)
        a = pdf.epw - 2 * pdf.c_margin - 1
        pdf.multi_cell(a, 10,
                       txt=section_content,
                       align='L', new_x=XPos.LEFT)
        pdf.image(image_name, w=150)


# Instantiation of inherited class
pdf = PDF()
pdf.alias_nb_pages()

# add outlines
pdf.set_font("Helvetica")
pdf.set_section_title_styles(
    # Level 0 titles:
    TitleStyle(
        font_family="Times",
        font_style="B",
        font_size_pt=24,
        color=128,
        underline=True,
        t_margin=10,
        l_margin=10,
        b_margin=0,
    ),
    # Level 1 subtitles:
    TitleStyle(
        font_family="Times",
        font_style="B",
        font_size_pt=20,
        color=128,
        underline=True,
        t_margin=10,
        l_margin=20,
        b_margin=5,
    ),
)
pdf.add_page()
pdf.set_y(50)
pdf.set_font(size=40)
new_section(pdf, "Doc Title", align="C")
pdf.set_font(size=12)
pdf.insert_toc_placeholder(render_toc)

# define sections and content
section_and_content = {'Vector model overview':
                           {'Vector model overview':
                                ['The EMOD vector model inherits the generic model functionality and introduces vector '
                                 'transmission and mosquito population dynamics. Interventions can be deployed within '
                                 'simulations for a variety of transmission settings with different transmission '
                                 'intensities, vector behaviors, and seasonally-driven ecologies. Climate data is '
                                 'necessary to simulate the effect of climatalogical impacts on vector biology. To use '
                                 'the vector model, set the configuration parameter Simulation_Type to VECTOR_SIM.\nThe'
                                 ' figure below demonstrates the main components of the vector EMOD simulation type.',
                                 'malariaSIR.png'
                                 ],
                            'Model implementation structure':
                                ['There are two categories of possible implementations of the basic model, each with different '
                                 'computational efficiencies, resolutions, and flexibilities. The first is an individual model, '
                                 'where it simulates every individual mosquito in the population or can utilize a sampled subset'
                                 ' of mosquitoes to represent the population as the whole. The second is a modified cohort '
                                 'simulation, with or without explicit mosquito ages.',
                                 None
                                 ]
                            },
                       'Malaria model':
                           {'Malaria model':
                                ['The malaria model inherits the functionality of the vector model and introduces human immunity,'
                                 ' within-host parasite dynamics, effects of antimalarial drugs, and other aspects of malaria '
                                 'biology to simulate malaria transmission. For example, individuals can have multiple infections'
                                 ' and both innate and adaptive responses to antigens. To use the malaria model, set the '
                                 'configuration parameter Simulation_Type to MALARIA_SIM.',
                                 None
                                 ],
                            'Model components':
                                ['The malaria model is complex, with numerous configurable parameters. '
                                 'The following network diagram breaks down the model into various model components, '
                                 'and illustrates how they interact with one another. The components on the network '
                                 'diagram correspond to the structural components listed below. Note that there is not '
                                 'perfect overlap between the labels on the network diagram and the structural '
                                 'components; this is because the network is drawn with increased detail in order to '
                                 'provide clarity in how the model functions and the components interact. The following '
                                 'pages will describe in detail how the structural components function.',
                                 'malaria_network_schematic.png'
                                 ]
                            }
                       }

for i, (section_tile, section_content) in enumerate(section_and_content.items()):
    pdf.start_section(f"{i + 1}. {section_tile}")
    subsection_names, text_and_images = section_content.keys(), section_content.values()
    for j, subsection_name, text_and_images in zip(range(len(subsection_names)), subsection_names, text_and_images):
        section_text, image_name = text_and_images
        pdf.start_section(f'{i + 1}.{j + 1} ' + subsection_name, level=1)
        new_section(
            pdf,
            section_text,
            image_name=image_name
        )
        pdf.add_page()

# pdf.add_page()
# for section_name, section_content, image_name in zip(section_names, section_content_list, image_names):
#     pdf.new_section(section_name=section_name,
#                     section_content=section_content,
#                     image_name=image_name)

pdf.output(pdf_name)
