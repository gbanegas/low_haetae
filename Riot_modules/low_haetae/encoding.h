#ifndef HAETAE_ENCODING_H
#define HAETAE_ENCODING_H

#include "params.h"
#include <stdint.h>
#include <stddef.h>

#define encode_h HAETAE_NAMESPACE(encode_h)
uint16_t encode_h(uint8_t *buf, const int32_t *h);
#define decode_h HAETAE_NAMESPACE(decode_h)
uint16_t decode_h(int32_t *h, const uint8_t *buf, uint16_t size_in);

/*
 * Low-SRAM helper: decode h to uint16_t (all HAETAE parameter sets keep h
 * in a small non-negative range).
 */
#define decode_h_u16 HAETAE_NAMESPACE(decode_h_u16)
uint16_t decode_h_u16(uint16_t *h, const uint8_t *buf, uint16_t size_in);
#define encode_hb_z1 HAETAE_NAMESPACE(encode_hb_z1)
uint16_t encode_hb_z1(uint8_t *buf, const int32_t *hb_z1);
#define decode_hb_z1 HAETAE_NAMESPACE(decode_hb_z1)
uint16_t decode_hb_z1(int32_t *hb_z1, const uint8_t *buf, uint16_t size_in);

/*
 * Low-SRAM helper: decode HB(z1) directly to int8_t symbols.
 * (HB(z1) is small; for HAETAE-120 this saves ~3KB vs int32_t output.)
 */
#define decode_hb_z1_i8 HAETAE_NAMESPACE(decode_hb_z1_i8)
uint16_t decode_hb_z1_i8(int8_t *hb_z1, const uint8_t *buf, uint16_t size_in);

#endif // HAETAE_ENCODING_H
