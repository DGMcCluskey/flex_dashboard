import tempfile
import matplotlib.pyplot as plt
from shiny import App, ui, render
import scanpy as sc

# Load your data
adata = sc.read_h5ad("data/flex_adata_clustered.h5ad")
genes = list(adata.var_names)
metadata = list(adata.obs.columns)

# UI
app_ui = ui.page_fluid(
    ui.h2("scRNA-seq Explorer"),
    ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_selectize("umap_gene", "Select gene", genes, multiple=False),
                ui.input_selectize("umap_meta", "Group by metadata", metadata, multiple=False),
            ),
            ui.output_plot("umap_meta"),
            ui.output_plot("umap_gene"),
        )
    ),
    ui.card(
        ui.layout_sidebar(
            ui.sidebar(
            ui.input_selectize("dotplot_gene", "Select gene", genes, multiple=False),
            ui.input_selectize("dotplot_meta", "Group by metadata", metadata, multiple=False)
        ),
        ui.output_image("dotplot"),  # dotplot displayed as image
        #ui.output_plot("violin"),
    )
)
)

# Server
def server(input, output, session):

    @output
    @render.plot
    def umap_meta():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_meta(), ax=ax, show=False)
        return fig

    @output
    @render.plot
    def umap_gene():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_gene(), ax=ax, show=False)
        return fig

    @output
    @render.image
    def dotplot():
        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)

    # 1. Create your own Matplotlib figure and axis
        fig, ax = plt.subplots(figsize=(6, 4))
    
    # 2. Draw the Scanpy dotplot on that axis
        sc.pl.dotplot(
            adata,
            var_names=[input.dotplot_gene()],
            groupby=input.dotplot_meta(),
            ax=ax,
            show=False  # needed to draw without popping up
        )
    # 3. Save the figure to disk
        fig.savefig(tmp.name, bbox_inches="tight")
        plt.close(fig)

    # 4. Return dict with "src" because Shiny requires this format
        return {"src": tmp.name}

# Run app
app = App(app_ui, server)