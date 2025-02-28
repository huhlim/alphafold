#!/usr/bin/env python
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Full AlphaFold protein structure prediction script."""
import json
import os
import pathlib
import pickle
import random
import sys
import time
from typing import Dict, Union, Optional, List

import jax
from absl import app
from absl import flags
from absl import logging
import numpy as np

from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import pipeline
from alphafold.data import pipeline_multimer
from alphafold.data import templates
from alphafold.data.tools import hhsearch
from alphafold.data.tools import hmmsearch
from alphafold.model import config
from alphafold.model import model
from alphafold.model import data
from alphafold.relax import relax

# Internal import (7716).

import libconfig_af
from libaf import *

# the main input arguments
flags.DEFINE_string(
    "fasta_path",
    None,
    "Paths to a FASTA file, If the FASTA file contains "
    "multiple sequences, then it will be folded as a multimer. ",
)

flags.DEFINE_string("data_dir", libconfig_af.data_dir, "Path to directory of supporting data.")
flags.DEFINE_string("output_dir", os.getcwd(), "Path to a directory that will store the results.")

# paths to executables
flags.DEFINE_string(
    "jackhmmer_binary_path", libconfig_af.jackhmmer_binary_path, "Path to the JackHMMER executable."
)
flags.DEFINE_string(
    "hhblits_binary_path", libconfig_af.hhblits_binary_path, "Path to the HHblits executable."
)
flags.DEFINE_string(
    "hhsearch_binary_path", libconfig_af.hhsearch_binary_path, "Path to the HHsearch executable."
)
flags.DEFINE_string(
    "hmmsearch_binary_path", libconfig_af.hmmsearch_binary_path, "Path to the hmmsearch executable."
)
flags.DEFINE_string(
    "hmmbuild_binary_path", libconfig_af.hmmbuild_binary_path, "Path to the hmmbuild executable."
)
flags.DEFINE_string(
    "kalign_binary_path", libconfig_af.kalign_binary_path, "Path to the Kalign executable."
)

# paths to databases
flags.DEFINE_string(
    "uniref90_database_path",
    libconfig_af.uniref90_database_path,
    "Path to the Uniref90 database for use by JackHMMER.",
)
flags.DEFINE_string(
    "mgnify_database_path",
    libconfig_af.mgnify_database_path,
    "Path to the MGnify database for use by JackHMMER.",
)
flags.DEFINE_string(
    "bfd_database_path",
    libconfig_af.bfd_database_path,
    "Path to the BFD database for use by HHblits.",
)
flags.DEFINE_string(
    "small_bfd_database_path",
    libconfig_af.small_bfd_database_path,
    'Path to the small version of BFD used with the "reduced_dbs" preset.',
)
flags.DEFINE_string(
    "uniref30_database_path",
    libconfig_af.uniref30_database_path,
    "Path to the uniref30 database for use by HHblits.",
)
flags.DEFINE_string(
    "uniprot_database_path",
    libconfig_af.uniprot_database_path,
    "Path to the Uniprot database for use by JackHMMer.",
)
flags.DEFINE_string(
    "pdb70_database_path",
    libconfig_af.pdb70_database_path,
    "Path to the PDB70 database for use by HHsearch.",
)
flags.DEFINE_string(
    "pdb_seqres_database_path",
    libconfig_af.pdb_seqres_database_path,
    "Path to the PDB seqres database for use by hmmsearch.",
)
flags.DEFINE_string(
    "template_mmcif_dir",
    libconfig_af.template_mmcif_dir,
    "Path to a directory with template mmCIF structures, each named <pdb_id>.cif",
)
flags.DEFINE_string(
    "max_template_date",
    libconfig_af.max_template_date,
    "Maximum template release date to consider. " "Important if folding historical test sets.",
)
flags.DEFINE_string(
    "obsolete_pdbs_path",
    libconfig_af.obsolete_pdbs_path,
    "Path to file containing a mapping from obsolete PDB IDs to the PDB IDs"
    "of their replacements.",
)

# presets
flags.DEFINE_enum(
    "db_preset",
    "full_dbs",
    ["full_dbs", "reduced_dbs"],
    "Choose preset MSA database configuration - "
    "smaller genetic database config (reduced_dbs) or "
    "full genetic database config  (full_dbs)",
)
flags.DEFINE_enum(
    "model_preset",
    "monomer",
    ["monomer", "monomer_casp14", "monomer_ptm", "multimer"],
    "Choose preset model configuration - the monomer model, "
    "the monomer model with extra ensembling, monomer model with "
    "pTM head, or multimer model",
)
flags.DEFINE_integer(
    "random_seed",
    None,
    "The random seed for the data "
    "pipeline. By default, this is randomly generated. Note "
    "that even if this is set, Alphafold may still not be "
    "deterministic, because processes like GPU inference are "
    "nondeterministic.",
)
flags.DEFINE_integer(
    "num_predictions_per_model",
    1,
    "How many "
    "predictions (each with a different random seed) will be "
    "generated per model. E.g. if this is 2 and there are 5 "
    "models then there will be 10 predictions per input. "
    "Note: this FLAG only applies if model_preset=multimer",
)
flags.DEFINE_boolean(
    "use_precomputed_msas",
    True,
    "Whether to read MSAs that "
    "have been written to disk. WARNING: This will not check "
    "if the sequence, database or configuration have changed.",
)

# custom arguments
flags.DEFINE_integer("cpu", 8, "Number of processors for sequence searches")
flags.DEFINE_boolean("jit", True, "compile using jax.jit")
flags.DEFINE_float(
    "max_sequence_identity", -1.0, "Maximum sequence identity for template prefilter"
)
flags.DEFINE_boolean("use_relax", True, "Whether to use AMBER local energy minimization")
flags.DEFINE_boolean("use_templates", True, "Whether to use PDB database")
flags.DEFINE_boolean("use_msa", True, "Whether to use MSA")
flags.DEFINE_list("template_mask", None, "Template masks. (e.g., 1-10,20-30)")
flags.DEFINE_boolean(
    "remove_msa_for_template_aligned", False, "Remove MSA information for template aligned region"
)
flags.DEFINE_integer("max_msa_clusters", None, "Number of maximum MSA clusters")
flags.DEFINE_integer("max_extra_msa", None, "Number of extra sequences")
flags.DEFINE_list("model_names", None, "Model configs to be run")
flags.DEFINE_list("msa_path", None, "User input MSA")
flags.DEFINE_list("pdb_path", None, "User input structure")
flags.DEFINE_boolean("unk_pdb", False, "Make input PDB residue names UNK")
flags.DEFINE_string("custom_templates", None, "User input templates")
flags.DEFINE_integer("num_recycle", 3, "The number of recycling")
flags.DEFINE_boolean("multimer", False, "Whether to use the multimer modeling hack")
flags.DEFINE_boolean("feature_only", False, "Whether to generate features.pkl only")
flags.DEFINE_boolean(
    "use_gpu_relax",
    False,
    "Whether to relax on GPU."
    "Relax on GPU can be much faster than CPU, so it is "
    "recommended to enable if possible. GPUs must be available"
    " if this setting is enabled.",
)

# presets
flags.DEFINE_boolean(
    "preset_refine",
    False,
    "Running refinement mode."
    'same as "--use_templates=False --use_msa=False --model_names=0 --unk_pdb=True'
    "requires --pdb_path",
)

FLAGS = flags.FLAGS

MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3


def _check_flag(flag_name: str, other_flag_name: str, should_be_set: bool):
    if should_be_set != bool(FLAGS[flag_name].value):
        verb = "be" if should_be_set else "not be"
        raise ValueError(
            f"{flag_name} must {verb} set when running with "
            f'"--{other_flag_name}={FLAGS[other_flag_name].value}".'
        )


def predict_structure(
    fasta_path: str,
    fasta_name: str,
    msa_path: Union[str, List[str]],
    pdb_path: Union[str, List[str]],
    output_dir_base: str,
    data_pipeline: Union[pipeline.DataPipeline, pipeline_multimer.DataPipeline],
    model_runners: Dict[str, model.RunModel],
    amber_relaxer: relax.AmberRelaxation,
    remove_msa_for_template_aligned: bool,
    template_mask: List[str],
    feature_only: bool,
    random_seed: int,
):
    """Predicts structure using AlphaFold for the given sequence."""

    logging.info("Predicting %s", fasta_name)
    timings = {}
    output_dir = os.path.join(output_dir_base, fasta_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    msa_output_dir = os.path.join(output_dir, "msas")
    if not os.path.exists(msa_output_dir):
        os.makedirs(msa_output_dir)

    # Get features.
    # modified to re-use features.pkl file, if it exists.
    t_0 = time.time()
    features_output_path = os.path.join(output_dir, "features.pkl")
    if os.path.exists(features_output_path):
        with open(features_output_path, "rb") as f:
            feature_dict = pickle.load(f)
    else:
        feature_dict = data_pipeline.process(
            input_fasta_path=fasta_path,
            input_msa_path=msa_path,
            input_pdb_path=pdb_path,
            msa_output_dir=msa_output_dir,
        )

        # Write out features as a pickled dictionary.
        with open(features_output_path, "wb") as f:
            pickle.dump(feature_dict, f, protocol=4)

    # apply the "remove_msa_for_template_aligned_regions" protocol
    if template_mask is not None:
        feature_dict = apply_template_mask(feature_dict, template_mask)
    if remove_msa_for_template_aligned:
        feature_dict = remove_msa_for_template_aligned_regions(feature_dict)

    timings["features"] = time.time() - t_0
    if feature_only:
        return

    unrelaxed_pdbs = {}
    relaxed_pdbs = {}
    ranking_confidences = {}

    # Run the models.
    num_models = len(model_runners)
    for model_index, (model_name, model_runner) in enumerate(model_runners.items()):
        unrelaxed_pdb_path = os.path.join(output_dir, f"unrelaxed_{model_name}.pdb")
        relaxed_output_path = os.path.join(output_dir, f"relaxed_{model_name}.pdb")
        if amber_relaxer:
            final_output_path = relaxed_output_path
        else:
            final_output_path = unrelaxed_pdb_path
        result_output_path = os.path.join(output_dir, f"result_{model_name}.pkl")
        if os.path.exists(final_output_path) and os.path.exists(result_output_path):
            # skip running this model and re-use pre-existing results.
            logging.info("Skipping %s", relaxed_output_path)
            with open(result_output_path, "rb") as fp:
                prediction_result = pickle.load(fp)
                ranking_confidences[model_name] = prediction_result["ranking_confidence"]

            with open(final_output_path) as fp:
                pdb_str = fp.read()
            if amber_relaxer:
                relaxed_pdbs[model_name] = pdb_str
            else:
                unrelaxed_pdbs[model_name] = pdb_str
            continue
        #
        logging.info("Running model %s on %s", model_name, fasta_name)

        t_0 = time.time()
        model_random_seed = model_index + random_seed * num_models
        processed_feature_dict = model_runner.process_features(
            feature_dict, random_seed=model_random_seed
        )
        timings[f"process_features_{model_name}"] = time.time() - t_0
        # processed_feat_path = os.path.join(output_dir, f"features_{model_name}.pkl")
        # with open(processed_feat_path, 'wb') as f:
        #    pickle.dump(processed_feature_dict, f, protocol=4)

        t_0 = time.time()
        prediction_result = model_runner.predict(
            processed_feature_dict, random_seed=model_random_seed
        )
        t_diff = time.time() - t_0
        timings[f"predict_benchmark_{model_name}"] = t_diff
        logging.info("Total JAX model %s on %s predict time: %.1fs", model_name, fasta_name, t_diff)

        plddt = prediction_result["plddt"]
        ranking_confidences[model_name] = prediction_result["ranking_confidence"]

        # Save the model outputs.
        with open(result_output_path, "wb") as f:
            pickle.dump(prediction_result, f, protocol=4)

        # retrieve custom features for outputs
        retrieve_custom_features(processed_feature_dict, feature_dict)

        # Add the predicted LDDT in the b-factor column.
        # Note that higher predicted LDDT value means higher model confidence.
        plddt_b_factors = np.repeat(plddt[:, None], residue_constants.atom_type_num, axis=-1)
        unrelaxed_protein = protein.from_prediction(
            features=processed_feature_dict,
            result=prediction_result,
            b_factors=plddt_b_factors,
            remove_leading_feature_dimension=not model_runner.multimer_mode,
        )

        unrelaxed_pdbs[model_name] = protein.to_pdb(unrelaxed_protein)
        with open(unrelaxed_pdb_path, "w") as f:
            f.write(unrelaxed_pdbs[model_name])

        # Relax the prediction.
        if amber_relaxer:
            t_0 = time.time()
            relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            timings[f"relax_{model_name}"] = time.time() - t_0
            relaxed_pdbs[model_name] = relaxed_pdb_str

            # Save the relaxed PDB.
            with open(relaxed_output_path, "w") as f:
                f.write(relaxed_pdb_str)

    # Rank by model confidence and write out relaxed PDBs in rank order.
    ranked_order = []
    for idx, (model_name, _) in enumerate(
        sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)
    ):
        ranked_order.append(model_name)
        ranked_output_path = os.path.join(output_dir, f"ranked_{idx}.pdb")
        with open(ranked_output_path, "w") as f:
            if amber_relaxer:
                f.write(relaxed_pdbs[model_name])
            else:
                f.write(unrelaxed_pdbs[model_name])

    ranking_output_path = os.path.join(output_dir, "ranking_debug.json")
    with open(ranking_output_path, "w") as f:
        label = "iptm+ptm" if "iptm" in prediction_result else "plddts"
        f.write(json.dumps({label: ranking_confidences, "order": ranked_order}, indent=4))

    logging.info("Final timings for %s: %s", fasta_name, timings)

    timings_output_path = os.path.join(output_dir, "timings.json")
    with open(timings_output_path, "w") as f:
        f.write(json.dumps(timings, indent=4))


def main(argv):
    if len(argv) > 1:
        raise app.UsageError("Too many command-line arguments.")
    #
    # modified: disabling JIT compilation
    if not FLAGS.jit:
        jax.config.update("jax_disable_jit", True)
    #
    # CHECK databases and executables
    for tool_name in ("jackhmmer", "hhblits", "hhsearch", "hmmsearch", "hmmbuild", "kalign"):
        if not FLAGS[f"{tool_name}_binary_path"].value:
            raise ValueError(
                f'Could not find path to the "{tool_name}" binary. '
                "Make sure it is installed on your system."
            )
    #
    if FLAGS.preset_refine:
        FLAGS.use_templates = False
        FLAGS.use_msa = False
        FLAGS.model_names = [0]
        FLAGS.unk_pdb = True
        if FLAGS.pdb_path is None:
            raise ValueError("--pdb_path is required for the refine preset")
    #
    use_small_bfd = FLAGS.db_preset == "reduced_dbs"
    if use_small_bfd:
        FLAGS.bfd_database_path = None
        FLAGS.uniref30_database_path = None
    else:
        FLAGS.small_bfd_database_path = None
    _check_flag("small_bfd_database_path", "db_preset", should_be_set=use_small_bfd)
    _check_flag("bfd_database_path", "db_preset", should_be_set=not use_small_bfd)
    _check_flag("uniref30_database_path", "db_preset", should_be_set=not use_small_bfd)

    run_multimer_system = "multimer" in FLAGS.model_preset
    if run_multimer_system:
        FLAGS.pdb70_database_path = None
    else:
        FLAGS.pdb_seqres_database_path = None
        FLAGS.uniprot_database_path = None
    _check_flag("pdb70_database_path", "model_preset", should_be_set=not run_multimer_system)
    _check_flag("pdb_seqres_database_path", "model_preset", should_be_set=run_multimer_system)
    _check_flag("uniprot_database_path", "model_preset", should_be_set=run_multimer_system)

    if FLAGS.msa_path:
        msa_path = []
        for pth in FLAGS.msa_path:
            if os.path.exists(pth):
                msa_path.append(pth)
            else:
                msa_path.append(None)
    else:
        msa_path = None
    if not run_multimer_system and (msa_path is not None):
        msa_path = msa_path[0]

    if FLAGS.pdb_path:
        pdb_path = []
        for pth in FLAGS.pdb_path:
            if os.path.exists(pth):
                pdb_path.append(pth)
            else:
                pdb_path.append(None)
    else:
        pdb_path = None
    if not run_multimer_system and (pdb_path is not None):
        pdb_path = pdb_path[0]

    if FLAGS.custom_templates is not None:
        FLAGS.template_mmcif_dir = FLAGS.custom_templates
        FLAGS.pdb70_database_path = "%s/pdb70" % FLAGS.custom_templates

    if FLAGS.multimer:
        FLAGS.use_templates = False
        if FLAGS.use_msa and msa_path is None:
            raise ValueError("The Multimer modeling hack requires an MSA input")
        if run_multimer_system:
            raise ValueError(
                "The Multimer modeling hack cannot be run with --model_preset=multimer"
            )

    if FLAGS.model_preset == "monomer_casp14":
        num_ensemble = 8
    else:
        num_ensemble = 1

    # PREPARE for running prediction
    fasta_name = pathlib.Path(FLAGS.fasta_path).stem

    # TEMPLATEs
    if FLAGS.use_templates:
        if run_multimer_system:
            template_searcher = hmmsearch.Hmmsearch(
                binary_path=FLAGS.hmmsearch_binary_path,
                hmmbuild_binary_path=FLAGS.hmmbuild_binary_path,
                database_path=FLAGS.pdb_seqres_database_path,
            )
            template_featurizer = templates.HmmsearchHitFeaturizer(
                mmcif_dir=FLAGS.template_mmcif_dir,
                max_template_date=FLAGS.max_template_date,
                max_hits=MAX_TEMPLATE_HITS,
                kalign_binary_path=FLAGS.kalign_binary_path,
                release_dates_path=None,
                obsolete_pdbs_path=FLAGS.obsolete_pdbs_path,
                max_sequence_identity=FLAGS.max_sequence_identity,
            )
        else:
            template_searcher = hhsearch.HHSearch(
                binary_path=FLAGS.hhsearch_binary_path, databases=[FLAGS.pdb70_database_path]
            )
            template_featurizer = templates.HhsearchHitFeaturizer(
                mmcif_dir=FLAGS.template_mmcif_dir,
                max_template_date=FLAGS.max_template_date,
                max_hits=MAX_TEMPLATE_HITS,
                kalign_binary_path=FLAGS.kalign_binary_path,
                release_dates_path=None,
                obsolete_pdbs_path=FLAGS.obsolete_pdbs_path,
                max_sequence_identity=FLAGS.max_sequence_identity,
            )
    else:
        template_searcher = None
        template_featurizer = None

    # Input Conformation
    if pdb_path is None:
        conformation_info_extractor = None
    else:
        conformation_info_extractor = templates.ConformationInfoExactractor(
            kalign_binary_path=FLAGS.kalign_binary_path, unk_pdb=FLAGS.unk_pdb
        )

    # PIPELINE
    monomer_data_pipeline = pipeline.DataPipeline(
        jackhmmer_binary_path=FLAGS.jackhmmer_binary_path,
        hhblits_binary_path=FLAGS.hhblits_binary_path,
        uniref90_database_path=FLAGS.uniref90_database_path,
        mgnify_database_path=FLAGS.mgnify_database_path,
        bfd_database_path=FLAGS.bfd_database_path,
        uniref30_database_path=FLAGS.uniref30_database_path,
        small_bfd_database_path=FLAGS.small_bfd_database_path,
        template_searcher=template_searcher,
        template_featurizer=template_featurizer,
        template_conformation=conformation_info_extractor,
        use_small_bfd=use_small_bfd,
        use_msa=FLAGS.use_msa,
        use_precomputed_msas=FLAGS.use_precomputed_msas,
        is_multimer=FLAGS.multimer,
        n_cpu=FLAGS.cpu,
    )
    num_predictions_per_model = FLAGS.num_predictions_per_model
    if run_multimer_system:
        data_pipeline = pipeline_multimer.DataPipeline(
            monomer_data_pipeline=monomer_data_pipeline,
            jackhmmer_binary_path=FLAGS.jackhmmer_binary_path,
            uniprot_database_path=FLAGS.uniprot_database_path,
            use_precomputed_msas=FLAGS.use_precomputed_msas,
            n_cpu=FLAGS.cpu,
        )
    else:
        data_pipeline = monomer_data_pipeline
    #
    if FLAGS.model_names is None:
        FLAGS.model_names = [0, 1, 2, 3, 4]
    else:
        FLAGS.model_names = [int(x) for x in FLAGS.model_names]

    model_runners = {}
    model_names = config.MODEL_PRESETS[FLAGS.model_preset]
    for i, model_name in enumerate(model_names):
        if i not in FLAGS.model_names:
            continue
        model_config = config.model_config(model_name)
        if run_multimer_system:
            model_config.model.num_ensemble_eval = num_ensemble
            model_config.model.num_recycle = FLAGS.num_recycle
        else:
            model_config.data.eval.num_ensemble = num_ensemble
            model_config.data.common.num_recycle = FLAGS.num_recycle
        #
        # modify MSA
        if FLAGS.max_msa_clusters is not None:
            model_config.data.eval.max_msa_clusters = FLAGS.max_msa_clusters
        if FLAGS.max_extra_msa is not None:
            model_config.data.common.max_extra_msa = FLAGS.max_extra_msa
        #
        model_params = data.get_model_haiku_params(model_name=model_name, data_dir=FLAGS.data_dir)
        model_runner = model.RunModel(model_config, model_params, jit_compile=FLAGS.jit)
        for i in range(num_predictions_per_model):
            model_runners[f"{model_name}_pred_{i}"] = model_runner

    logging.info("Have %d models: %s", len(model_runners), list(model_runners.keys()))

    # RELAX
    if FLAGS.use_relax:
        amber_relaxer = relax.AmberRelaxation(
            max_iterations=RELAX_MAX_ITERATIONS,
            tolerance=RELAX_ENERGY_TOLERANCE,
            stiffness=RELAX_STIFFNESS,
            exclude_residues=RELAX_EXCLUDE_RESIDUES,
            max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
            use_gpu=FLAGS.use_gpu_relax,
        )
    else:
        amber_relaxer = None

    # RANDOMSEED
    random_seed = FLAGS.random_seed
    if random_seed is None:
        random_seed = random.randrange((sys.maxsize) // len(model_runners))
    logging.info("Using random seed %d for the data pipeline", random_seed)

    # RUN PREDICTION
    predict_structure(
        fasta_path=FLAGS.fasta_path,
        fasta_name=fasta_name,
        msa_path=msa_path,
        pdb_path=pdb_path,
        output_dir_base=FLAGS.output_dir,
        data_pipeline=data_pipeline,
        model_runners=model_runners,
        amber_relaxer=amber_relaxer,
        remove_msa_for_template_aligned=FLAGS.remove_msa_for_template_aligned,
        template_mask=FLAGS.template_mask,
        feature_only=FLAGS.feature_only,
        random_seed=random_seed,
    )


if __name__ == "__main__":
    flags.mark_flags_as_required(
        [
            "fasta_path",
        ]
    )

    app.run(main)
