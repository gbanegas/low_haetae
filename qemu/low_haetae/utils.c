#include "utils.h"


/*************************************************
 * Name:        generate_seed_from_one_source
 *
 * Description: Use SHAKE256 to generate a seed from one byte source.
 *
 * Arguments:   - uint8_t *seed:    pointer to buffer of seed to be generated
 *              - size_t seed_len:  desired length of the seed
 *              - uint8_t *src:    pointer to source to be absorbed
 *              - size_t src_len:  length of source
 *
 * Returns void
 **************************************************/
void generate_seed_from_one_source(uint8_t *seed, size_t seed_len,
		const uint8_t *src, size_t src_len) {
	xof256_state state;
	xof256_absorbe_once(&state, src, src_len);
	xof256_squeeze(seed, seed_len, &state);
	return;
}