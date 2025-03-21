import multiprocessing

import pickle
import tqdm
import pandas as pd
import numpy as np
import argparse

# import time
# from datetime import datetime

import os

try:
    from pyensembl import EnsemblRelease
except ImportError:
    print("Use genotype_selection env, you need pyensembl with GRCh37 available")
    exit(0)


def map_range(x, data, ensmbl_field="Field", chromosome_dict={}):
    """
    Maps a string ENSEMBL Range (ENSBL1-ENSBL2) to a set of ENSBML Terms.
    A dictionary called "chromosome_dict" must be defined in enclosing namespace!
    :param x: Pandas dataset row. String column 'gid' in format "{ENSBL1}-{ENSBL2}" must exist.
    :param data: PyEnsembl Datasource
    :return: String encoding of set of strings of ENSBML terms in range @x['gid']

    originally by Simon Sasse, slightly adapted

    """
    # in case we dont have a range, just return
    field_range = x[0]  # just for the error message

    try:
        if not ("-" in x[ensmbl_field]):
            # print(f'single ensid: {x[ensmbl_field]}')
            return x[ensmbl_field]
        field_range = x[ensmbl_field]
        chrom = x["CHROM"]  # Extract chromosome number
        # print("Chrom: ", chrom)
        genes = field_range.split("-")
        # print("Genes ", genes)
        start_pos = 0 if genes[0] == "CHR_START" else data.gene_by_id(genes[0]).start
        end_pos = 1000000000 if genes[1] == "CHR_END" else data.gene_by_id(genes[1]).end
        if (
            not chrom in chromosome_dict
        ):  # If the current chromosome data hasn't been loaded
            # print("inside IF")
            all_g = data.genes(contig=int(chrom))
            all_g_pos = [g.start for g in all_g]
            all_g_end = [g.end for g in all_g]
            all_g_ensg = [g.gene_id for g in all_g]
            df = pd.DataFrame(
                data={"pos": all_g_pos, "end": all_g_end, "gname": all_g_ensg}
            )
            chromosome_dict[chrom] = df
        in_range = chromosome_dict[chrom][
            np.logical_and(
                chromosome_dict[chrom].pos >= start_pos,
                chromosome_dict[chrom].pos <= end_pos,
            )
        ]
        return str(set([in_range.loc[e].gname for e in in_range.index]))

    except (IndexError, TypeError) as inst:
        print("Error with:", field_range)
        print(inst)
        return str(set())


def conversion(x):
    """
    Function that every process executes. Take in a str and process it to a numpy array and saves it to memorymaps
    :param x: (tuple) nr of position and lines provided from file
    :return: mean
    """
    # TODO: Replacing N/A with 0 dosage, don't forget this
    _, info = x
    files, rownrs, out = info
    array = np.asarray(
        [
            item.replace("NA", "0")
            for sublist in [liste.split("\t")[6:] for liste in out]
            for item in sublist
        ],
        dtype=np.float16,
    )

    if args.compress:
        array = np.ceil(array * 127).astype("uint8")

        for nr, file in enumerate(files):
            readmap = np.memmap(
                f"{file}.X.npy",
                dtype="uint8",
                mode="r+",
                order="C",
                shape=(array.shape[0]),
                offset=array.shape[0] * rownrs[nr],
            )
            readmap[:] = array
            # reshape? mean pro spalte
            del readmap

    else:
        for nr, file in enumerate(files):
            readmap = np.memmap(
                f"{file}.X.npy",
                dtype="float16",
                mode="r+",
                order="C",
                shape=(array.shape[0]),
                offset=array.shape[0] * 2 * 1 * rownrs[nr],
            )

            readmap[:] = array
            # reshape? mean pro spalte
            del readmap
    return array.reshape(-1, length).mean(axis=0)


def read_chroms_expanded(out_path, file_list, selections):
    """
    Generator that returns file lines
    :param out_path: outputpath of memmaps
    :param file_list: files to read
    :param selections:
    :return: memmapfiles, row to write to and current line
    """
    relevant_rows = set()
    for sel in selections:
        relevant_rows = set(sel.ravel()).union(relevant_rows)

    # write dict that tracks which line goes to which file
    # write row tracking dict that tracks how many lines have been written to a file
    file_dict = {name: [] for name in relevant_rows}
    row_tracker = {
        f'{out_path}/splits/all_chrs_{(".".join(file.split(".")[:-1])).split("/")[-1].replace("controlled", args.suffix)}': 0
        for file in sel_dict.keys()
    }

    for file, array in sel_dict.items():
        for name in array:
            name = int(name)
            if name in file_dict.keys():
                file_name = (".".join(file.split(".")[:-1])).split("/")[-1]
                file_name = file_name.replace("controlled", args.suffix)
                file_name = f"all_chrs_{file_name}"
                file_dict[name].append(f"{out_path}/splits/{file_name}")

    # open all files and skip header
    files_open = [open(file, "r") for file in file_list]
    header = [file.readline() for file in files_open]
    nr = 0
    while True:
        lines = [file.readline() for file in files_open]

        if nr in file_dict.keys():
            outfiles = file_dict[nr]

            row_nrs = []
            for file in outfiles:
                row_nrs.append(row_tracker[file])
                row_tracker[file] += 1

            yield outfiles, row_nrs, lines

        nr += 1

        if "" in lines:
            files_closed = [file.close() for file in files_open]
            break


def to_npy_split(
    var_path,
    raw_path,
    out_path,
    pickpath,
    eid_mapping=None,
    eid_from_col=None,
    eid_to_col=None,
):
    if eid_mapping is not None:
        assert eid_from_col is not None
        assert eid_to_col is not None
        eid_mapping = (
            pd.read_csv(eid_mapping, sep="\t")
            .set_index(eid_from_col)[eid_to_col]
            .to_dict()
        )

    if out_path.endswith(".X.npy"):
        out_path = out_path[:-6]
    os.makedirs(out_path + "/splits", exist_ok=True)

    splits_path = out_path + "/splits"

    write_any_x = False
    for file_path in args.pickpath:
        file_name = (".".join(file_path.split(".")[:-1])).split("/")[-1]
        file_name = file_name.replace("controlled", args.suffix)
        file_name = f"all_chrs_{file_name}"
        if not (os.path.isfile(f"{out_path}/splits/{file_name}.X.npy")):
            # and os.path.isfile(f'{out_path}/splits/{file_name}.obs.csv')):
            write_any_x = True
            print(f"did NOT find {out_path}/splits/{file_name}.X.npy")
        else:
            print(f"found {out_path}/splits/{file_name}.X.npy")

    if not write_any_x:
        print("All .X.npy files present already, not writing them again.")

    # for var

    # determine the variants present in the raw data after removing duplicates by plink
    variants = None

    valid_raw_files = []
    for chrom in range(1, 23):
        try:
            tmp = pd.read_csv(raw_path.format(chrom), sep="\t", nrows=1)
        except FileNotFoundError:
            print("Didn't find " + raw_path.format(chrom))
            continue
        except pd.errors.EmptyDataError:
            print("Empty file: " + raw_path.format(chrom))
            continue
        valid_raw_files.append(raw_path.format(chrom))
        try:
            variants = variants + list(tmp.columns[6:])
        except TypeError:
            variants = list(tmp.columns[6:])

    del tmp
    global length
    length = len(variants)
    print(f"Found {len(variants)} variants in .raw files.")

    # for obs (not optimized)
    raw_df = pd.read_csv(valid_raw_files[0], sep="\t", usecols=list(range(5)))
    if eid_mapping is not None:
        obs = pd.DataFrame(
            index=raw_df["IID"].map(eid_mapping)
        )  # .loc[ex_raw["IID"] > 0, "IID"])
    else:
        obs = pd.DataFrame(index=raw_df["IID"])  # .loc[ex_raw["IID"] > 0, "IID"])
    obs.index.name = "index"
    if not np.issubdtype(obs.index.dtype, np.number) and "_" in obs.index[0]:
        obs.index = obs.index.map(lambda x: x.split("_")[0])
        print(
            "Found underscores in IID column, removed everything after the underscore"
        )

    # for file_path in args.pickpath:
    #     file_name = (".".join(file_path.split(".")[:-1])).split('/')[-1]
    #     file_name = file_name.replace('controlled', args.suffix)
    #     file_name = f'all_chrs_{file_name}'
    #
    #     obs_opath = f'{out_path}/splits/{file_name}.obs.csv'
    obs_opath = f"{out_path}/obs.csv"
    obs.to_csv(obs_opath)
    print(f"Wrote obs to {obs_opath} with {len(obs)} samples")
    del raw_df

    selections = []
    global sel_dict
    sel_dict = {}
    sel_dict_raw = {}
    obs.index = obs.index.map(str)

    # read the pick_list files for selection. safe info in different formats
    for file in pickpath:

        with open(file, "r") as f:
            pick_list = f.read().splitlines()

        # selraw.append(pick_list)
        selections.append(np.argwhere(obs.index.isin(pick_list)))
        sel_dict[file] = np.argwhere(obs.index.isin(pick_list))
        sel_dict_raw[file] = obs.index.isin(pick_list)

        print(f"Read pick file {file} with {len(pick_list)} entries.")
        # print(pick_list[:5])t
        # print(str(sel_dict[file])[:500])

    if write_any_x:
        print(f"Starting to write X.")

        #  prep memmaps
        for name in sel_dict.keys():
            if sel_dict[name].shape[0] == 0 or len(variants) == 0:
                raise RuntimeError(
                    "No variants or samples in selection, probably a mapping error for eids or variants, check obs and var files"
                )
            file_name = (".".join(name.split(".")[:-1])).split("/")[-1]
            file_name = file_name.replace("controlled", args.suffix)
            file_name = f"all_chrs_{file_name}"
            print(
                f"Creating memmap for {name} with shape {(sel_dict[name].shape[0], len(variants))} at: {out_path}/splits/{file_name}.X.npy"
            )
            if args.compress:
                out_mmapped = np.memmap(
                    f"{out_path}/splits/{file_name}.X.npy",
                    dtype="uint8",
                    mode="w+",
                    order="C",
                    shape=(sel_dict[name].shape[0], len(variants)),
                )
                out_mmapped[:] = 0
                del out_mmapped
            else:
                out_mmapped = np.memmap(
                    f"{out_path}/splits/{file_name}.X.npy",
                    dtype="float16",
                    mode="w+",
                    order="C",
                    shape=(sel_dict[name].shape[0], len(variants)),
                )
                out_mmapped[:] = 0
                del out_mmapped

        it = read_chroms_expanded(out_path, valid_raw_files, selections)  # Generator

        # Multiprocessing starts here:
        print("Starting the multiprocessing")
        with multiprocessing.Pool() as p:

            # looks weird and is weird. but enables tqdm for nice overview of progress.
            # provides lines to the conversion processes
            mean = np.zeros(length)
            c = 0
            for x in tqdm.tqdm(p.imap_unordered(conversion, enumerate(it))):
                mean += x
                c += 1
            mean = mean / c
            # print(mean,mean.shape)
            # print(c)

        print(f"Done writing X to {x_opath} with shape {(len(obs), len(variants))}")
    else:
        mean = np.nan
    # read the actual variant annotation
    if args.skipvar:
        print("Not making a var.csv file")
        var_opath = out_path + ".var.csv"
        var = pd.read_csv(var_path)

    else:
        pre_var = pd.read_csv(var_path)
        assert (
            pre_var["ID"].duplicated().sum() == 0
        ), "Found duplicated IDs in var file, results will be misleading"

        pre_var = pre_var.sort_values(by=["CHROM", "POS"])

        pre_var["_raw_id1"] = pre_var.apply(
            axis=1, func=lambda x: f'{x["ID"]}_{x["ALT"]}'
        )
        pre_var["_raw_id2"] = pre_var.apply(
            axis=1,
            func=lambda x: f'{x["CHROM"]}:{x["POS"]}_{x["REF"]}_{x["ALT"]}_{x["ALT"]}',
        )
        pre_var["_raw_id3"] = pre_var.apply(
            axis=1, func=lambda x: f'{x["ID"]}_{x["REF"]}(/{x["ALT"]})'
        )

        pre_var = pre_var[
            np.logical_or(
                pre_var["_raw_id1"].isin(variants),
                np.logical_or(
                    pre_var["_raw_id2"].isin(variants),
                    pre_var["_raw_id3"].isin(variants),
                ),
            )
        ]

        # TODO: Verify order?

        # put down the ensemblids
        chromosome_dict = {}
        datasource = EnsemblRelease(75)  # this is the version we need for GRCh37.5

        pre_var["Field"] = pre_var.apply(
            axis=1,
            func=lambda x: map_range(
                x,
                datasource,
                ensmbl_field="ANN_Gene_ID",
                chromosome_dict=chromosome_dict,
            ),
        )
        pre_var["ValueType"] = "Genetic"
        pre_var["WasProcessed"] = "False"
        pre_var[
            "GeneticFrequency"
        ] = mean  # TODO: Danger! mean is currently only calculated for the in the selection included id

        pre_var = pre_var.rename(
            columns={
                "CHROM": "GeneticChromosome",
                "POS": "GeneticPosition",
                "REF": "GeneticRef",
                "ALT": "GeneticAlt",
                "ANN_Annotation_Impact": "GeneticImpact",
                "ANN_Gene_Name": "GeneticGName",
                "ANN_Transcript_BioType": "GeneticTranscript",
                "ANN_Annotation": "GeneticEffect",
            }
        )

        var = pre_var
        var_opath = out_path + ".var.csv"
        var.to_csv(var_opath)
        print(f"Wrote var to {var_opath} with {len(var)} variants")

    for name in sel_dict.keys():

        file_name = (".".join(name.split(".")[:-1])).split("/")[-1]
        file_name = file_name.replace("controlled", args.suffix)
        file_name = f"all_chrs_{file_name}"

        # this is just a copy operation

        var_opath = os.path.join(splits_path, file_name + ".var.csv")
        # print(var_opath)
        var.to_csv(var_opath)
        print(f"Wrote var to {var_opath} with {len(var)} variants")

        # define the selection index

        obs.index = obs.index.map(str)
        selection = sel_dict_raw[name]

        obs_opath = os.path.join(splits_path, file_name + ".obs.csv")
        obs[selection].to_csv(obs_opath)
        print(f"Wrote obs to {obs_opath} with {len(obs)} samples")

        # uns.pkl file

        uns = {}
        if args.compress:
            uns["dtype"] = "uint8_127"
        else:
            uns["dtype"] = "float16"

        uns_opath = os.path.join(splits_path, file_name + ".uns.pkl")
        with open(uns_opath, "wb") as f:
            pickle.dump(uns, f)
        print(f'Wrote uns to {uns_opath} with dtype {uns["dtype"]} ')


if __name__ == "__main__":
    file_path_parser = argparse.ArgumentParser(description="adjust_annotated_file")
    # arguments
    file_path_parser.add_argument("--i_var", help="Input path to future var dataframe.")
    file_path_parser.add_argument(
        "--i_raw", help="Input path to raw genotype data from plink."
    )
    file_path_parser.add_argument("--output", help="Output path.")
    file_path_parser.add_argument("--suffix", help="suffix for file naming")
    file_path_parser.add_argument(
        "--pickpath", help="path to IDs used for splitting.", nargs="+"
    )
    file_path_parser.add_argument(
        "--compress", help="compression to uint8", action="store_true"
    )
    file_path_parser.add_argument(
        "--skipvar", help="skip writing of var file", action="store_true"
    )
    file_path_parser.add_argument(
        "--eid_mapping_file", help="table mapping eids", required=False
    )
    file_path_parser.add_argument("--eid_map_from", help="Column name", required=False)
    file_path_parser.add_argument("--eid_map_to", help="Column name", required=False)

    args = file_path_parser.parse_args()

    x_opath = args.output + ".X.npy"

    to_npy_split(
        args.i_var,
        args.i_raw.replace("chr_1", "chr_{}"),
        args.output,
        args.pickpath,
        args.eid_mapping_file,
        args.eid_map_from,
        args.eid_map_to,
    )  # TODO: this might not  be ideal when some other path is specified
    print("Done")