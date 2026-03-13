import shinyswatch
import tempfile
import matplotlib.pyplot as plt
from shiny import App, ui, render
import scanpy as sc
import numpy as np
import marsilea as ma
import marsilea.plotter as mp
import io

# Load your data
adata = sc.read_h5ad("data/flex_adata_clustered.h5ad") 
genes = list(adata.var_names)
metadata = list(adata.obs.columns)
metadata = [col for col in metadata if col in ["sample_id", "pool_id", "poms_id", 
"celltypist_label_foetal_lung_celltypist", "leiden_res_1"]]

# UI
app_ui = ui.page_fluid(
    ui.h2("FLEX data dashboard"),
    #row 1
    ui.layout_columns(
        ui.card(
        ui.card_header("UMAP by metadata"),
        ui.output_plot("umap_meta"),
        ui.input_radio_buttons("umap_meta_input", "Group by metadata", 
        choices=metadata, selected="sample_id"),
        ui.download_button("download_umap_meta", "Save as PDF")
    ),

    ui.card(
        ui.card_header("Overlay gene expression"),
        ui.output_plot("umap_gene"),
        ui.input_selectize("umap_gene_input", "Select gene", genes, selected="PTPRC"),
        ui.output_text_verbatim("gene_summary"),
        ui.download_button("download_umap_gene_expression", "Save as PDF")
    ),
    col_widths=(8,4),
),
    #row 2
    ui.layout_columns(
    ui.card(
        ui.card_header("Detailed gene expression"),
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_radio_buttons("dotplot_meta_input", "Group by metadata", 
                choices =metadata, selected="sample_id"),
                ui.input_selectize("dotplot_gene_input", "Select gene", genes, multiple=True, 
                selected=["PTPRC", "CD3D", "LYZ"]),
                ui.download_button("download_dotplot", "Save as PDF")
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

        genes = list(input.dotplot_gene_input())
        groupby = input.dotplot_meta_input()

        if len(genes) == 0:
            return None

        tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)

        # -----------------------------
        # Aggregate expression
        # -----------------------------

        agg = sc.get.aggregate(
            adata[:, genes],
            by=groupby,
            func=["mean", "count_nonzero"]
        )
        sc.pp.scale(agg, layer="mean")  

        # number of cells per group
        agg.obs["cell_counts"] = adata.obs[groupby].value_counts()

        # -----------------------------
        # Extract matrices for Marsilea
        # -----------------------------

        mean_expr = agg.layers["mean"]
        frac_expr = agg.layers["count_nonzero"] / agg.obs["cell_counts"].values[:, None]

        # -----------------------------
        # Build Marsilea dotplot
        # -----------------------------

        s = ma.SizedHeatmap(
            size=frac_expr,
            color=mean_expr,
            cmap="Reds",
            edgecolor="black",
            width=8,
            height=5,
            size_legend_kws=dict(
            colors="#538bbf",
            title="Fraction of cells\nin groups (%)",
            labels=["20%", "40%", "60%", "80%", "100%"],
            show_at=[0.2, 0.4, 0.6, 0.8, 1.0],
        ),
        color_legend_kws=dict(title="Mean expression\nin group"),
    )

        s.add_left(mp.Labels(agg.obs_names))
        s.add_bottom(mp.Labels(agg.var_names))
        s.add_legends()
        s.add_left(mp.Numbers(agg.obs["cell_counts"], color="#EEB76B", label="Count"))

        s.render()

        fig = plt.gcf()
        fig.savefig(tmp.name, bbox_inches="tight")
        plt.close(fig)

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

    @output
    @render.download()
    def download_umap_meta():
        buf = io.BytesIO()
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_meta_input(), ax=ax, show=False)
        fig.savefig(buf, format="pdf", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf, "umap_meta.pdf"

    @output
    @render.download()
    def download_umap_gene_expression():
        buf = io.BytesIO()
        fig, ax = plt.subplots()
        sc.pl.umap(adata, color=input.umap_meta_input(), ax=ax, show=False)
        fig.savefig(buf, format="pdf", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf, "umap_gene.pdf"

    @output
    @render.download()
    def download_dotplot():
        genes = list(input.dotplot_gene_input())
        if len(genes) == 0:
            return None

        fig, ax = plt.subplots(figsize=(10,5))
        
        # Marsilea or Scanpy plotting code
        # e.g. s.render() or sc.pl.dotplot(...)
        
        fig.savefig(buf := io.BytesIO(), format="pdf", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf, "dotplot.pdf"

# Run app
app = App(app_ui, server)