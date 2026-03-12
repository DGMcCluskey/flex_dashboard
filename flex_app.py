import shinyswatch
import tempfile
import matplotlib.pyplot as plt
from shiny import App, ui, render
import scanpy as sc
import numpy as np

# Load your data
adata = sc.read_h5ad("data/flex_adata_clustered.h5ad") 
genes = list(adata.var_names)
metadata = list(adata.obs.columns)
metadata = [col for col in metadata if col in ["sample_id", "pool_id", "poms_id", 
"celltypist_label_foetal_lung_celltypist", "leiden_res_1"]]
# UI
app_ui = ui.page_fluid(
    ui.h2("scRNA-seq Explorer"),
    #row 1
    ui.layout_columns(
        ui.card(
        ui.card_header("UMAP by metadata"),
        ui.output_plot("umap_meta"),
        ui.input_checkbox_group("umap_meta_input", "Group by metadata", 
        choices=metadata, selected="sample_id")
    ),

    ui.card(
        ui.card_header("Overlay gene expression"),
        ui.output_plot("umap_gene"),
        ui.input_selectize("umap_gene_input", "Select gene", genes, selected="PTPRC"),
        ui.output_text_verbatim("gene_summary")
    ),
    col_widths=(8,4),
),
    #row 2
    ui.layout_columns(
    ui.card(
        ui.card_header("Detailed gene expression"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_checkbox_group("dotplot_meta_input", "Group by metadata", 
                choices =metadata, selected="sample_id"),
                ui.input_selectize("dotplot_gene_input", "Select gene", genes, multiple=True, 
                selected=["PTPRC", "CD3D", "LYZ"])
            ),
            ui.output_plot("dotplot"),
        )
    )
    ),
    theme=shinyswatch.theme.solar
)

# Server
def server(input, output, session):

    @output
    @render.plot
    def umap_meta():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_meta_input(), ax=ax, show=False, 
        legend_fontsize=6, size =5)
        return fig

    @output
    @render.plot
    def umap_gene():
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_gene_input(), ax=ax, show=False, cmap = "magma_r", vmax="p99",
        size=5)
        return fig

    @output
    @render.image
    def dotplot():

        genes = list(input.dotplot_gene_input())  # convert tuple → list

    # if nothing selected, don't try to plot
        if len(genes) == 0:
            return None

        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)

        fig, ax = plt.subplots(figsize=(10,5))

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

    @output
    @render.text
    def gene_summary():

        gene = input.umap_gene_input()

        if gene is None:
            return ""

        x = adata[:, gene].X

    # handle sparse matrices (very common in scRNA)
        if hasattr(x, "toarray"):
            x = x.toarray()

        x = x.flatten()

        summary = f"""
    Cells expressing: {(x>0).sum()}
    Mean expression: {x.mean():.3f}
    Median expression: {float(np.median(x)):.3f}
    Max expression: {x.max():.3f}
    """

        return summary

# Run app
app = App(app_ui, server)