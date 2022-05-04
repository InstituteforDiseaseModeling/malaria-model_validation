from fpdf import FPDF

title = 'Malaria Model Overview'

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
        self.set_font('Arial', 'B', 25)
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
        self.cell(w, 9, title, 1, 1, 'C', 1)
        # Line break
        self.ln(20)

    # Page footer
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')

# Instantiation of inherited class
pdf = PDF()
pdf.alias_nb_pages()
pdf.add_page()
pdf.set_font('Times', '', 24)
pdf.cell(0, 10, 'Model components', 0, 1)

pdf.set_font('Times', '', 12)
pdf.multi_cell(0, 10, 'The malaria model is complex, with numerous configurable parameters. The following network diagram breaks down the model into various model components, and illustrates how they interact with one another. The components on the network diagram correspond to the structural components listed below. Note that there is not perfect overlap between the labels on the network diagram and the structural components; this is because the network is drawn with increased detail in order to provide clarity in how the model functions and the components interact. The following pages will describe in detail how the structural components function.', 0, 1)
pdf.image('malaria_network_schematic.png', w=150)


pdf.set_font('Times', '', 24)
pdf.cell(0, 10, 'Vector model overview', 0, 1)

pdf.set_font('Times', '', 12)
pdf.multi_cell(0, 10, 'The EMOD vector model inherits the generic model functionality and introduces vector transmission and mosquito population dynamics. Interventions can be deployed within simulations for a variety of transmission settings with different transmission intensities, vector behaviors, and seasonally-driven ecologies. Climate data is necessary to simulate the effect of climatalogical impacts on vector biology. To use the vector model, set the configuration parameter Simulation_Type to VECTOR_SIM.', 0, 1)
pdf.multi_cell(0, 10, 'The figure below demonstrates the main components of the vector EMOD simulation type.', 0, 1)
pdf.image('malariaSIR.png', w=150)

pdf.output('Malaria_model_validation_output.pdf', 'F')