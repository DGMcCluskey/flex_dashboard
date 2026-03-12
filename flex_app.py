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
    #row 1
    ui.layout_columns(
        
        ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize("umap_meta_input", "Group by metadata", metadata, multiple=False),
                ),
                ui.output_plot("umap_meta"),
            )
        ),
        ui.card(
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_selectize("umap_gene_input", "Select gene", genes, multiple=False) 
                ),
                ui.output_plot("umap_gene")
            )
        )
    ),
    #row 2
    ui.layout_columns(
    ui.card(
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_selectize("dotplot_meta_input", "Group by metadata", metadata, multiple=False),
                ui.input_selectize("dotplot_gene_input", "Select gene", genes, multiple=True)
            ),
            ui.output_plot("dotplot"),
        )
    )
    )
)

# Server
def server(input, output, session):

    @output
    @render.plot
    def umap_meta():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_meta_input(), ax=ax, show=False)
        return fig

    @output
    @render.plot
    def umap_gene():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_gene_input(), ax=ax, show=False)
        return fig

    @output
    @render.image
    def dotplot():

        genes = list(input.dotplot_gene_input())  # convert tuple → list

    # if nothing selected, don't try to plot
        if len(genes) == 0:
            return None

        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)

        fig, ax = plt.subplots(figsize=(6, 4))

        sc.pl.dotplot(
            adata,
            var_names=genes,
            groupby=input.dotplot_meta_input(),
            ax=ax,
            show=False
        )

        fig.savefig(tmp.name, bbox_inches="tight")
        plt.close(fig)

    # 4. Return dict with "src" because Shiny requires this format
        return {"src": tmp.name}

# Run app
app = App(app_ui, server)