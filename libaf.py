#!/usr/bin/env python

def remove_msa_for_template_aligned_regions(feature_dict):
    if 'template_all_atom_masks' in feature_dict:
        mask = feature_dict['template_all_atom_masks']
    elif 'template_all_atom_mask' in feature_dict:
        mask = feature_dict['template_all_atom_mask']
    mask = (mask.sum(axis=(0,2)) > 0)
    #
    # need to check further for multimer_mode
    if 'deletion_matrix_int' in feature_dict:
        feature_dict['deletion_matrix_int'][:,mask] = 0
    else:
        feature_dict['deletion_matrix'][:,mask] = 0
    feature_dict['msa'][:,mask] = 21
    return feature_dict

def retrieve_custom_features(processed_feature_dict, feature_dict):
    for name in ['for_pdb_record']:
        if name in feature_dict:
            processed_feature_dict[name] = feature_dict[name]

