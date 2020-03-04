//
// Created by Andrew Bailey on 2019-07-19.
//

#include "VariantCall.hpp"

/**
 * Constructor for VariantCall.
 * Sets contig, strand, ref index and bases
 */
VariantCall::VariantCall(string contig1, string strand1, uint64_t reference_index1, string bases1) :
    contig(std::move(contig1)), strand(std::move(strand1)), reference_index(reference_index1), bases(std::move(bases1)) {}
/**
 * Constructor for VariantCall.
 * Sets reference index
 */
VariantCall::VariantCall() : reference_index(-1) {}

/**
 * Deconstructor for VariantCall. Need to close file
 */
VariantCall::~VariantCall()
= default;


