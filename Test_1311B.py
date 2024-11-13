from shiny import App, render, ui, reactive
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

# Load the Excel file and the required sheets
file_path_heatmap = "heatmaps_recall50.xlsx"
xls = pd.ExcelFile(file_path_heatmap)

# Load pathway names and group classifications
pathway_names_df = xls.parse("pathway_names_classes")
# Rename columns for clarity
pathway_names_df.columns = ["pathway_id", "pathway_name", "unknown", "group"]

# Clean pathway names to remove any special characters and extra spaces
pathway_names_df["pathway_name"] = pathway_names_df["pathway_name"].str.replace("\xa0", " ", regex=True).str.strip()
pathway_names_df["group"] = pathway_names_df["group"].str.replace("\xa0", " ", regex=True).str.strip()

# Define the pathway groups for filtering
pathway_classes = [
    "Membrane transport",
    "Signal transduction",
    "Folding, sorting and degradation",
    "Translation",
    "Carbohydrate metabolism",
    "Energy metabolism",
    "Glycan biosynthesis and metabolism",
    "Metabolism of other amino acids"
]

# Default drugs for Cell Wall and Protein Synthesis, which users can override
default_cell_wall_drugs = ["BAC", "CSD", "VAN"]
default_protein_synthesis_drugs = ["CHL", "DOX", "ERY", "TET"]

# Define available sheets for methods
sheet_options = ["chemogenomics", "ml", "docking", "transcriptomics", "metabolomics"]

# Define custom color mapping based on the legend from the reference image
def custom_colormap(value):
    if value == -100:
        return "#FFC0CB"  # Light Red
    elif value == -200:
        return "#FFA07A"  # Light Blue
    elif value == -300:
        return "#FF6347"  # Medium Blue
    elif value == -400:
        return "#FF4500"  # Darker Blue
    elif -1000 <= value <= -500:
        return "#B22222"  # Black for 21-50 hits
    elif -2000 <= value <= -1100:
        return "#8B0000"
    elif value == -1:
        return "#E0FFE0"  # Light Blue
    elif value == -2:
        return "#ADFF2F"  # Medium Blue
    elif value == -3:
        return "#7FFF00"  # Darker Blue
    elif value == -4:
        return "#32CD32"  # Deep Blue
    elif -10 <= value <= -5:
        return "#228B22"
    elif -20 <= value <= -11:
        return "#006400"
    elif value == 1:
        return "#E0FFFF"  # Light Blue
    elif value == 2:
        return "#ADD8E6"
    elif value == 3:
        return "#87CEEB"
    elif value == 4:
        return "#4169E1"
    elif 5 <= value <= 10:
        return "#0000CD"
    elif 11 <= value <= 50:
        return "#000000"
    elif value == 0:
        return "white"
    return "gray"

# Define display transformation for negative values
def transform_display_value(value):
    if value < 0:
        # For red legend (e.g., -100, -200, etc.), display as 1, 2, etc.
        if value % 100 == 0:
            return abs(value // 100)
        # For green legend (e.g., -1, -2, etc.), display as 1, 2, etc.
        else:
            return abs(value)
    return value  # Return as-is for non-negative values

# Define UI for the application with multi-select inputs for drugs
app_ui = ui.page_fluid(
    ui.div(
        ui.h2("Ultra-high-throughput screening of antimicrobial combination therapies using a two-stage transparent machine learning model"),
        ui.h4("Margaret M. Reuter, Katherine Lev, Jon Albo, Harkirat Singh Arora, Nemo Liu, Madeline Shay, Debmalya Sarkar, Aaron Robida, David H. Sherman, Rudy J. Richardson, Nate J. Cira, Sriram Chandrasekaran"),
        style="text-align: center;"
    ),
    ui.div(
        ui.a("Visit Systems Biology Lab", href="https://systemsbiologylab.org/", target="_blank", style="position: absolute; top: 10px; right: 10px; color: blue; font-size: 16px;"),
    ),
    ui.input_select("method", "Choose Method(s):", choices=sheet_options, multiple=True),
    ui.input_select("pathway_class", "Choose Pathway Class(es):", choices=pathway_classes, multiple=True),
    
    # Multi-select dropdown for Cell Wall Drugs
    ui.input_select("cell_wall_drugs", "Choose Cell Wall Drugs:", choices=default_cell_wall_drugs, multiple=True),
    
    # Multi-select dropdown for Protein Synthesis Drugs
    ui.input_select("protein_synthesis_drugs", "Choose Protein Synthesis Drugs:", choices=default_protein_synthesis_drugs, multiple=True),
    
    ui.output_plot("heatmap_plot", height="800px")
)

# Define server logic
def server(input, output, session):
    # Reactive expression to load data based on selected method, pathway class, and drugs
    @reactive.Calc
    def selected_data():
        methods = input.method()
        pathway_classes_selected = input.pathway_class()

        # Get selected drugs dynamically
        selected_cell_wall_drugs = input.cell_wall_drugs() if input.cell_wall_drugs() else default_cell_wall_drugs
        selected_protein_synthesis_drugs = input.protein_synthesis_drugs() if input.protein_synthesis_drugs() else default_protein_synthesis_drugs

        # Filter pathways based on selected classes
        selected_pathways = pathway_names_df[pathway_names_df["group"].isin(pathway_classes_selected)]["pathway_name"].tolist()

        # Load and merge data from selected methods
        combined_data = []
        for method in methods:
            data = xls.parse(method)
            data.columns = ["pathway_name"] + list(data.columns[1:])
            data["pathway_name"] = data["pathway_name"].str.replace("\xa0", " ", regex=True).str.strip()

            # Filter data based on selected pathways
            data = data[data["pathway_name"].isin(selected_pathways)]
            combined_data.append(data)
        
        # If no data, return None
        if not combined_data:
            print("Selected data is empty. Please check the dataset or select a different method.")
            return None, None, None, None
        
        # Concatenate data from all selected methods
        data = pd.concat(combined_data, ignore_index=True)

        # Separate data by drug targets and select only the specified drugs
        cell_wall_data = data[data.columns.intersection(["pathway_name"] + list(selected_cell_wall_drugs))].set_index("pathway_name")
        protein_synthesis_data = data[data.columns.intersection(["pathway_name"] + list(selected_protein_synthesis_drugs))].set_index("pathway_name")
        
        # Fill any NaN values with zeros
        cell_wall_data = cell_wall_data.fillna(0)
        protein_synthesis_data = protein_synthesis_data.fillna(0)

        # Generate titles based on selected methods and pathways
        method_title = ", ".join(methods)
        pathway_title = ", ".join(pathway_classes_selected)
        
        return cell_wall_data, protein_synthesis_data, method_title, pathway_title, selected_cell_wall_drugs, selected_protein_synthesis_drugs

    # Render the composite heatmap plot
    @output
    @render.plot
    def heatmap_plot():
        cell_wall_data, protein_synthesis_data, method_title, pathway_title, selected_cell_wall_drugs, selected_protein_synthesis_drugs = selected_data()

        if cell_wall_data is None or protein_synthesis_data is None:
            # Return an empty plot if data is invalid
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'Error: Selected data is empty or not available.', ha='center', va='center')
            return fig

        # Load the legend image
        legend_img = mpimg.imread("C:/Users/debss/Documents/shiny/legend.png")

        # Create the figure with an additional subplot for the legend
        fig, axs = plt.subplots(1, 3, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 1, 0.5]})

        # Plot heatmaps on the first two subplots
        for ax, data, title_prefix, selected_drugs in zip(
            [axs[0], axs[1]], 
            [cell_wall_data, protein_synthesis_data], 
            ["Drug Target: Cell Wall", "Drug Target: Protein Synthesis"],
            [selected_cell_wall_drugs, selected_protein_synthesis_drugs]
        ):
            ax.set_title(f"{title_prefix} ({method_title}) - {pathway_title}")
            
            # Center x-ticks and y-ticks
            ax.set_xticks(np.arange(len(selected_drugs)) + 0.5)  # Shift x-tick positions to center
            ax.set_xticklabels(selected_drugs, ha='center', va='top')
            
            ax.set_yticks(np.arange(len(data.index)) + 0.5)  # Shift y-tick positions to center
            ax.set_yticklabels(data.index, ha='right', va='center')

            for (y, x), value in np.ndenumerate(data.values):
                color = custom_colormap(value)
                rect = plt.Rectangle((x, y), 1, 1, facecolor=color, edgecolor="black", linewidth=1)
                ax.add_patch(rect)
                
                # Display transformed value based on the logic
                display_value = transform_display_value(value)
                ax.text(x + 0.5, y + 0.5, f"{display_value}", ha='center', va='center', fontsize=9,
                        color="black" if color != "white" else "black")

            # Set grid limits and invert y-axis for orientation
            ax.set_xlim(0, len(selected_drugs))
            ax.set_ylim(0, len(data.index))
            ax.invert_yaxis()

        # Add the legend image on the third subplot
        axs[2].imshow(legend_img)
        axs[2].axis('off')  # Hide axes for the legend image

        # Adjust layout to prevent overlap
        fig.tight_layout()

        return fig

# Define the ASGI app object
app = App(app_ui, server)
