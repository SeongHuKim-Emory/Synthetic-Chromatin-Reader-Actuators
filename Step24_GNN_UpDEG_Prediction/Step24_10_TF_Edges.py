import polars as pl
from pathlib import Path
import sys
import polars.selectors as cs


def main():
    """
    Main function to process ChIP-seq data and augment the graph edges.
    """
    
    # --- 1. Import files ---
    
    # Define base paths
    try:
        BASE_DIR = Path(__file__).resolve().parent
    except NameError:
        # Fallback for interactive environments (like Jupyter)
        BASE_DIR = Path.cwd()


    PUBLIC_DIR = BASE_DIR / "../../public"
    MULTIOMICS_DIR = PUBLIC_DIR / "Multiomics/Step24_GNN_UpDEG_Prediction"
    CHIPSEQ_DIR = PUBLIC_DIR / "Public Dataset/ChIPseq"


    # Define file paths
    nodes_path = MULTIOMICS_DIR / "Nodes.csv"
    edges_path = MULTIOMICS_DIR / "Edges.csv"


    # List of all ChIP-seq BED files to process
    chipseq_filenames = [
        "MCF7 ARID3A ChIPseq GSE91595 ENCFF823HUJ hg38.bed",
        "MCF7 ATF7 ChIPseq GSE106025 ENCFF261DZR hg38.bed",
        "MCF7 BMI1 ChIPseq GSE105933 ENCFF846CPP hg38.bed",
        "MCF7 CEBPB ChIPseq GSM1010889 ENCFF240ZHW hg38.bed",
        "MCF7 CHD1 ChIPseq GSE91698 ENCFF587BZB hg38.bed",
        "MCF7 CLOCK ChIPseq GSE127438 ENCFF594YBN hg38.bed",
        "MCF7 COPS2 ChIPseq GSE105356 ENCFF395PJC hg38.bed",
        "MCF7 CREB1 ChIPseq GSE105525 ENCFF599FYT hg38.bed",
        "MCF7 CSDE1 ChIPseq GSE105291 ENCFF576UBB hg38.bed",
        "MCF7 CTBP1 ChIPseq GSE91938 ENCFF882ZZZ hg38.bed",
        "MCF7 CTCF ChIPseq GSE123219 ENCFF278FNP hg38.bed",
        "MCF7 CUX1 ChIPseq GSE91415 ENCFF184PHK hg38.bed",
        "MCF7 DDX20 ChIPseq GSE105517 ENCFF565XMS hg38.bed",
        "MCF7 DPF2 ChIPseq GSE91598 ENCFF870ORM hg38.bed",
        "MCF7 E2F8 ChIPseq GSE127615 ENCFF096JXF hg38.bed",
        "MCF7 E4F1 ChIPseq GSE127582 ENCFF434IDL hg38.bed",
        "MCF7 EGR1 ChIPseq GSM1010844 ENCFF427RSA hg38.bed",
        "MCF7 ELF1 ChIPseq GSE105503 ENCFF989DLF hg38.bed",
        "MCF7 ELK1 ChIPseq GSE91713 ENCFF427XNV hg38.bed",
        "MCF7 EP300 ChIPseq GSM1010800 ENCFF290DJX hg38.bed",
        "MCF7 ESRRA ChIPseq GSE92187 ENCFF214ICQ hg38.bed",
        "MCF7 FOS ChIPseq GSE105734 ENCFF519VZH hg38.bed",
        "MCF7 FOXA1 ChIPseq GSE105305 ENCFF112JVK hg38.bed",
        "MCF7 FOXK2 ChIPseq GSE91785 ENCFF222VDW hg38.bed",
        "MCF7 FOXM1 ChIPseq GSM1010769 ENCFF563MHZ hg38.bed",
        "MCF7 GABPA ChIPseq GSM1010864 ENCFF678DJM hg38.bed",
        "MCF7 GATA3 ChIPseq GSE127656 ENCFF217KJR hg38.bed",
        "MCF7 GATAD2B ChIPseq GSE91723 ENCFF401YUZ hg38.bed",
        "MCF7 GTF2F1 ChIPseq GSE91869 ENCFF819AEG hg38.bed",
        "MCF7 HCFC1 ChIPseq GSE91992 ENCFF796EMB hg38.bed",
        "MCF7 HDAC2 ChIPseq GSM1010825 ENCFF661VXV hg38.bed",
        "MCF7 HDGF ChIPseq GSE91567 ENCFF480XAU hg38.bed",
        "MCF7 HES1 ChIPseq GSE105287 ENCFF132MJX hg38.bed",
        "MCF7 HSF1 ChIPseq GSE91444 ENCFF766UTT hg38.bed",
        "MCF7 JUN ChIPseq GSE91550 ENCFF769TUY hg38.bed",
        "MCF7 JUND ChIPseq GSM1010892 ENCFF398OOD hg38.bed",
        "MCF7 LARP7 ChIPseq GSE105864 ENCFF734FBR hg38.bed",
        "MCF7 MAFK ChIPseq GSE127366 ENCFF425GWD hg38.bed",
        "MCF7 MAX ChIPseq GSM1010863 ENCFF624CRN hg38.bed",
        "MCF7 MAZ ChIPseq GSE91633 ENCFF290CNC hg38.bed",
        "MCF7 MBD2 ChIPseq GSE96480 ENCFF853ULZ hg38.bed",
        "MCF7 MLLT1 ChIPseq GSE91754 ENCFF972FNO hg38.bed",
        "MCF7 MNT ChIPseq GSE91968 ENCFF655RVZ hg38.bed",
        "MCF7 MTA1 ChIPseq GSE91687 ENCFF070JYN hg38.bed",
        "MCF7 MTA2 ChIPseq GSE91864 ENCFF556NWZ hg38.bed",
        "MCF7 MTA3 ChIPseq GSE91727 ENCFF975UOJ hg38.bed",
        "MCF7 NBN ChIPseq GSE91899 ENCFF065ASC hg38.bed",
        "MCF7 NCOA3 ChIPseq GSE105681 ENCFF778ZBT hg38.bed",
        "MCF7 NEUROD1 ChIPseq GSE127346 ENCFF796KBR hg38.bed",
        "MCF7 NFIB ChIPseq GSE105751 ENCFF079JBB hg38.bed",
        "MCF7 NFRKB ChIPseq GSE105388 ENCFF532IPJ hg38.bed",
        "MCF7 NFXL1 ChIPseq GSE105498 ENCFF142NAV hg38.bed",
        "MCF7 NONO ChIPseq GSE92159 ENCFF115WXJ hg38.bed",
        "MCF7 NR2F2 ChIPseq GSM1010837 ENCFF386FQQ hg38.bed",
        "MCF7 NRF1 ChIPseq GSE91522 ENCFF970MPY hg38.bed",
        "MCF7 PAX8 ChIPseq GSE105903 ENCFF490MXJ hg38.bed",
        "MCF7 PKNOX1 ChIPseq GSE92210 ENCFF043NUC hg38.bed",
        "MCF7 PML ChIPseq GSM1010838 ENCFF569JWM hg38.bed",
        "MCF7 POLR2A ChIPseq GSM822295 ENCFF235UTX hg38.bed",
        "MCF7 PPP1R10 ChIPseq GSE105558 ENCFF428GSG hg38.bed",
        "MCF7 RAD21 ChIPseq GSM1010791 ENCFF728MPA hg38.bed",
        "MCF7 RAD51 ChIPseq GSE105597 ENCFF062KLW hg38.bed",
        "MCF7 RCOR1 ChIPseq GSE91726 ENCFF160KBR hg38.bed",
        "MCF7 REST ChIPseq GSM1010891 ENCFF680JMZ hg38.bed",
        "MCF7 RFX1 ChIPseq GSE91448 ENCFF782YIG hg38.bed",
        "MCF7 RFX5 ChIPseq GSE105874 ENCFF753LKY hg38.bed",
        "MCF7 SIN3A ChIPseq GSE91789 ENCFF402DKZ hg38.bed",
        "MCF7 SIX4 ChIPseq GSE91630 ENCFF728MNA hg38.bed",
        "MCF7 SMARCA5 ChIPseq GSE105722 ENCFF526KAO hg38.bed",
        "MCF7 SMARCE1 ChIPseq GSE105546 ENCFF889QAN hg38.bed",
        "MCF7 SNIP1 ChIPseq GSE105188 ENCFF285KBA hg38.bed",
        "MCF7 SP1 ChIPseq GSE92014 ENCFF932UJX hg38.bed",
        "MCF7 SREBF1 ChIPseq GSE91561 ENCFF571VYR hg38.bed",
        "MCF7 SRF ChIPseq GSM1010839 ENCFF131TYZ hg38.bed",
        "MCF7 SUZ12 ChIPseq GSE105981 ENCFF488VNT hg38.bed",
        "MCF7 TAF1 ChIPseq GSM1010811 ENCFF075TIU hg38.bed",
        "MCF7 TARDBP ChIPseq GSE105812 ENCFF666QVW hg38.bed",
        "MCF7 TCF7L2 ChIPseq GSM816438 ENCFF332AZX hg38.bed",
        "MCF7 TCF12 ChIPseq GSM1010861 ENCFF987QLI hg38.bed",
        "MCF7 TEAD4 ChIPseq GSM1010860 ENCFF751VAZ hg38.bed",
        "MCF7 TRIM22 ChIPseq GSE127607 ENCFF238QAG hg38.bed",
        "MCF7 YBX1 ChIPseq GSE92114 ENCFF745UTS hg38.bed",
        "MCF7 ZBTB1 ChIPseq GSE105482 ENCFF192BPQ hg38.bed",
        "MCF7 ZBTB7B ChIPseq GSE105418 ENCFF339BDC hg38.bed",
        "MCF7 ZBTB11 ChIPseq GSE95952 ENCFF328ZJV hg38.bed",
        "MCF7 ZBTB40 ChIPseq GSE91661 ENCFF160TJG hg38.bed",
        "MCF7 ZFX ChIPseq GSE105562 ENCFF861DOL hg38.bed",
        "MCF7 ZHX2 ChIPseq GSE96441 ENCFF443WEJ hg38.bed",
        "MCF7 ZKSCAN1 ChIPseq GSE91769 ENCFF856VYW hg38.bed",
        "MCF7 ZNF8 ChIPseq GSE127545 ENCFF435UDQ hg38.bed",
        "MCF7 ZNF24 ChIPseq GSE105531 ENCFF258VTL hg38.bed",
        "MCF7 ZNF207 ChIPseq GSE91475 ENCFF761KBK hg38.bed",
        "MCF7 ZNF217 ChIPseq GSE105662 ENCFF197CVM hg38.bed",
        "MCF7 ZNF444 ChIPseq GSE127367 ENCFF896QHA hg38.bed",
        "MCF7 ZNF507 ChIPseq GSE127652 ENCFF987CPY hg38.bed",
        "MCF7 ZNF512B ChIPseq GSE105987 ENCFF320KEW hg38.bed",
        "MCF7 ZNF574 ChIPseq GSE127638 ENCFF281CNU hg38.bed",
        "MCF7 ZNF579 ChIPseq GSE105224 ENCFF223BRJ hg38.bed",
        "MCF7 ZNF592 ChIPseq GSE91425 ENCFF566KDV hg38.bed",
        "MCF7 ZNF687 ChIPseq GSE92150 ENCFF609HVM hg38.bed"
    ]


    print("Loading initial Nodes and Edges files...")
    # Load Nodes.csv
    try:
        nodes_df = pl.read_csv(nodes_path)
        nodes_df = nodes_df.with_columns(pl.col("node_id").cast(pl.Int64))
    except Exception as e:
        print(f"Error loading {nodes_path}: {e}")
        sys.exit(1)


    # Load Edges.csv
    try:
        edges_df = pl.read_csv(edges_path)
        edges_df = edges_df.with_columns(
            cs.integer().cast(pl.Int64)
        )
    except Exception as e:
        print(f"Error loading {edges_path}: {e}")
        sys.exit(1)


    print("Initial files loaded.")


    # --- 4. Prepare data for efficient lookups ---


    # Filter for promoter nodes
    promoter_types = [
        "promoter_proximal_downstream",
        "promoter_core",
        "promoter_proximal_upstream"
    ]
    promoter_nodes = nodes_df.filter(
        pl.col("node_type").is_in(promoter_types)
    ).select(["node_id", "chr", "start", "end"])


    # Create a dictionary for fast TF 'gene' node_id lookup
    gene_nodes = nodes_df.filter(pl.col("node_type") == "gene").select(["identifier", "node_id"])
    tf_to_node_id = {row[0]: row[1] for row in gene_nodes.rows()}


    # Get the set of original 'is_' columns from Edges.csv
    original_is_cols = {col for col in edges_df.columns if col.startswith("is_")}
    
    # --- 5, 6, 7. Process all ChIP-seq files ---
    
    all_new_edges_list = []
    all_new_tf_cols_set = set()
    
    current_max_edge_id = edges_df["edge_id"].max()
    new_edge_counter = 0


    print("Starting ChIP-seq file processing...")
    for filename in chipseq_filenames:
        try:
            tf_name = filename.split(" ")[1]
            new_col_name = f"is_{tf_name}"

            if tf_name not in tf_to_node_id:
                print(f"Warning: TF {tf_name} not found in Nodes file. Skipping.")
                continue
            
            all_new_tf_cols_set.add(new_col_name)
            
            tf_node_id = tf_to_node_id[tf_name]

            bed_path = CHIPSEQ_DIR / filename
            bed_df = pl.read_csv(
                bed_path,
                separator="\t",
                has_header=False,
                columns=[0, 1, 2],
                new_columns=["chr", "peak_start", "peak_end"],
                use_pyarrow=True,
                schema_overrides={"chr": pl.Utf8, "peak_start": pl.Int64, "peak_end": pl.Int64}
            )

            candidates = promoter_nodes.join(bed_df, on="chr", how="inner")
            
            found_overlaps = candidates.filter(
                (pl.col("start") < pl.col("peak_end")) & (pl.col("end") > pl.col("peak_start"))
            )

            overlapping_promoter_ids = found_overlaps.select("node_id").unique()

            num_new_edges = overlapping_promoter_ids.height
            if num_new_edges > 0:
                new_ids = pl.arange(
                    current_max_edge_id + 1 + new_edge_counter,
                    current_max_edge_id + 1 + new_edge_counter + num_new_edges,
                    eager=True
                )
                
                new_rows_df = pl.DataFrame({
                    "edge_id": new_ids,
                    "node_id_source": tf_node_id,
                    "node_id_target": overlapping_promoter_ids["node_id"],
                    new_col_name: 1
                })
                
                all_new_edges_list.append(new_rows_df)
                new_edge_counter += num_new_edges
            
            print(f"Processed {tf_name}: found {num_new_edges} new edges.")

        except Exception as e:
            print(f"Error processing file {filename}: {e}. Skipping.")

    print("All files processed. Combining DataFrames...")

    # --- 8. Combine original and new edges ---

    if not all_new_edges_list:
        print("No new edges were found. Exporting original Edges file.")
        final_edges_df = edges_df
    else:
        combined_new_edges = pl.concat(all_new_edges_list, how="diagonal")

        cols_to_add_to_new = [
            pl.lit(0, dtype=pl.Int64).alias(col) for col in original_is_cols 
            if col not in combined_new_edges.columns
        ]
        if cols_to_add_to_new:
            combined_new_edges = combined_new_edges.with_columns(cols_to_add_to_new)

        cols_to_add_to_old = [
            pl.lit(0, dtype=pl.Int64).alias(col) for col in all_new_tf_cols_set 
            if col not in original_is_cols
        ]
        if cols_to_add_to_old:
            edges_df = edges_df.with_columns(cols_to_add_to_old)
            
        combined_new_edges = combined_new_edges.fill_null(0)

        all_final_columns = edges_df.columns
        combined_new_edges = combined_new_edges.select(all_final_columns)

        # --- FIX: Cast all integer columns to Int64 ---
        # This ensures the schema matches the original edges_df before concatenation.
        combined_new_edges = combined_new_edges.with_columns(
            cs.integer().cast(pl.Int64)
        )

        final_edges_df = pl.concat([edges_df, combined_new_edges], how="vertical")
        print("DataFrames combined successfully.")

    # --- 9. Export the final Edges table ---
    output_path = MULTIOMICS_DIR / "Edges.csv"
    try:
        final_edges_df.write_csv(output_path)
        print(f"Successfully exported updated edges to {output_path}")
        print(f"Original edges: {edges_df.height}, New edges: {new_edge_counter}, Total edges: {final_edges_df.height}")
    except Exception as e:
        print(f"Error writing final CSV to {output_path}: {e}")


if __name__ == "__main__":
    main()


