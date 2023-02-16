# Running Baysor with type annotations

import Baysor
import Colors
import Plots
import CSV

using Distances
using ProgressMeter
using DataFrames
using DataFramesMeta
using NearestNeighbors
using RCall
using Statistics
using StatsBase

coords_pxl_csv = snakemake.input["coords_pxl_csv"]
centroid_profiles_csv = snakemake.input["centroid_profiles_csv"]
output_anno_csv = snakemake.output["anno_csv"]
output_time_json = snakemake.output["time_json"]

B = Baysor;

timings = Dict()

@info "Reading data"
# read rna
df_spatial, gene_names = B.load_df(coords_pxl_csv);

# read mean expression
centroid_df = CSV.read(centroid_profiles_csv, DataFrame);
centroid_exprs = Matrix(centroid_df[:, Symbol.(gene_names)]);
@info "Done reading data"

@info "Starting single molecule annotation"
# annotation transfer
start = time()
adjacent_points, adjacent_weights = B.build_molecule_graph(df_spatial, filter=false); # adjacency
B.append_confidence!(df_spatial, nn_id=16); # adding confidence

type_transfer_init = B.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    do_maximize=false, max_iters=1000, cell_type_exprs=centroid_exprs, n_iters_without_update=20);

timings["annotation"] = time()-start
@info "Done with single molecule annotation"
df_spatial[!, :clusterannotation] = centroid_df.ClusterName[type_transfer_init.assignment];

# Cell segmentation
## parameters come from the Baysor repository (default values or specialized values for the osm FISH dataset) (Baysor v0.5.2)
args = Dict([
    ("scale-std", "25%"), # default value from Baysor/src/cli/main.jl
    ("iters", 500), # from Baysor/examples/osm-FISH/README.md
    ("scale", 70.0), # from Baysor/configs/osm_fish.toml
    ("new-component-fraction", 0.3), # default value from Baysor/src/cli/common.jl
    ("num-cells-init", 30000), # from Baysor/examples/osm-FISH/README.md
    ("prior-segmentation-confidence", 0.2), # default value from Baysor/src/cli/main.jl
    ("min-molecules-per-cell", 30), # from Baysor/configs/osm_fish.toml
    ("new-component-weight", 0.2), # default value from Baysor/src/cli/common.jl
])
comp_genes = Vector{Int}[]
@info "Starting cell segmentation"
start = time()
bm_data = B.initialize_bmm_data(df_spatial; scale=args["scale"], scale_std=args["scale-std"],
        n_cells_init=args["num-cells-init"], prior_seg_confidence=args["prior-segmentation-confidence"],
        min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=0,
        adjacent_points=adjacent_points, adjacent_weights=adjacent_weights, na_genes=Vector{Int}(vcat(comp_genes...)));
history_depth = round(Int, args["iters"] * 0.1)
bm_data = B.bmm!(bm_data; n_iters=args["iters"], new_component_frac=args["new-component-fraction"], new_component_weight=args["new-component-weight"],
        min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth);

timings["segmentation"] = time()-start
@info "Done with cell segmentation"

# Save results
segmentated_df = B.get_segmentation_df(bm_data, convert(Vector{String},gene_names))
@info "Saving results to $(output_anno_csv)"
CSV.write(output_anno_csv, segmentated_df[!, ["cell","clusterannotation"]], compress=endswith(output_anno_csv,".gz"));

open(output_time_json, "w") do io
    write(io, "{\"baysor\":$(timings["annotation"]),\"baysor_segmenting_baysor\":$(timings["segmentation"])}")
end;

@info "All done!"
