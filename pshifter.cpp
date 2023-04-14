#include <complex>  // std::complex
#include <vector>   // std::vector
#include <cmath>
#include <cstdint>

#include <stdio.h>
#include <stdlib.h> // atof()
#include <cstring>

#define DELTA_T 0.05            // seconds
#define ROUND_2_INT(f) ((int)(f >= 0. ? (f + 0.5) : (f - 0.5)))

FILE *fp;
unsigned short CHANNELS;
unsigned int SAMPLE_RATE;       // Hz
unsigned short SAMPLE_SIZE;     // bytes
bool BYTE_FORMAT;               // 0 - signed, 1 - unsigned
bool BYTES_ORDER;               // 0 - little endian, 1 - big endian
unsigned int AUDIO_SIZE;        // bytes
unsigned short HDR_SIZE;        // size of audio file (bytes) except for AUDIO_SIZE

// =============================================== FFT BLOCK ===============================================

unsigned int findMSB(int x)
{
    unsigned int p = 0;

    while (x != 1) {
        x >>= 1;
        ++p;
    }

    return p;
}

int bitr(unsigned int x, unsigned int nb)
{
    x = ( x               << 16) | ( x               >> 16);
    x = ((x & 0x00FF00FF) <<  8) | ((x & 0xFF00FF00) >>  8);
    x = ((x & 0x0F0F0F0F) <<  4) | ((x & 0xF0F0F0F0) >>  4);
    x = ((x & 0x33333333) <<  2) | ((x & 0xCCCCCCCC) >>  2);
    x = ((x & 0x55555555) <<  1) | ((x & 0xAAAAAAAA) >>  1);
    // for 32 bit integer
    return ((x >> (32 - nb)) & (0xFFFFFFFF >> (32 - nb)));
}

std::vector<std::complex<float>> fft1d(const std::vector<std::complex<float>> &xi, char dir)
{
    unsigned int cnt = (unsigned int)xi.size();
    float nrm = 1.f / sqrt(cnt);
    unsigned int msb = findMSB(cnt);
    std::vector<std::complex<float>> xo(cnt);

    // pre-process the input data
    for (unsigned int j = 0; j != cnt; ++j)
        xo[j] = nrm * xi[bitr(j, msb)];

    cnt>>=1;
    int bm, bw;
    int i1, i2, j;
    std::complex<float> tmp;
    float ang, angc = dir*M_PI;

    // fft passes
    for (unsigned int i = 0; i != msb; ++i) {
        bm = 1 << i;    // butterfly mask
        bw = bm << 1;   // butterfly width
        ang = angc/bm;  // precomputation

        // fft butterflies
        for (j = 0; j != cnt; ++j) {
            i1 = ((j >> i) << (i + 1)) + (j & (bm-1));  // left wing.
            i2 = i1 ^ bm;                               // right wing
            tmp = xo[i1];

            xo[i1]+= xo[i2] * std::polar(1.f, ang * (i1 ^ bw));
            xo[i2] = xo[i2] * std::polar(1.f, ang * (i2 ^ bw)) + tmp;
        }
    }

    return xo;
}

// =================================================== PARSING HEADER ======================================================

unsigned int readHeader()
{
    unsigned int itype;
    unsigned short stype;
    unsigned short pos = 36;
    char ctype[4];

    if (fread(&itype, 1, 4, fp) != 4)                   { puts("Invalid file!\n"); return 0; }
    // RIFF or RIFX
    if (itype != 1179011410 && itype != 1481001298)     { puts("Invalid file format!\n"); return 0; }

    memcpy(ctype, &itype, 4);
    ctype[4] = 0;                       // delete newline character
    printf("File Type: %s\n", ctype);

    fread(&itype, 1, 4, fp);
    printf("File size: %u bytes\n", itype);

    fseek(fp, 20, 0);
    fread(&stype, 1, 2, fp);
    if (stype == 1 || stype == 65534) printf("Type of Format: (integer) PCM\n");
    else                            { printf("Invalid audio format!\n"); return 0; }

    fread(&stype, 1, 2, fp);
    printf("Channels: %u\n", stype);
    CHANNELS = stype;

    fread(&itype, 1, 4, fp);
    printf("Sample rate: %u\n", itype);
    SAMPLE_RATE = itype;

    fseek(fp, 34, 0);
    fread(&stype, 1, 2, fp);
    printf("Sample size: %u\n", stype);
    SAMPLE_SIZE = stype>>3;

    // 8 bit sample size is unsigned int, others are signed int
    if (stype == 8) {
        printf("Byte Format: Unsigned\n");
        BYTE_FORMAT = 1;
    }
    else {
        printf("Byte Format: Signed\n");
        BYTE_FORMAT = 0;
    }

    // RIFF header means little endian, RIFX header means big endian (?)
    if     (ctype[3] == 'F') {
        printf("Byte order: Little Endian\n");
        BYTES_ORDER = 0;
    }
    else /*(ctype[3] == 'X')*/ {
        printf("Byte order: Big Endian\n");
        BYTES_ORDER = 1;
    }

    // FINDING "data" SUBCHUNK
    itype = 0;
    do {
        if (fread(&itype, 1, 4, fp) != 4) { puts("Invalid file: no 'data' chunk found!\n"); return 0; }
        pos += 4;
    } while (itype != 1635017060);      // char[4] ("data") equals to (uint) 1635017060

    fread(&itype, 1, 4, fp);
    printf("Audio data size: %u bytes\n", itype);
    AUDIO_SIZE = itype;

    pos += 4;
    HDR_SIZE = pos;
    printf("Header size: %u bytes\n", HDR_SIZE);

    if (AUDIO_SIZE == 0) { puts("Invalid audio data size!\n"); return 0; }
    while (AUDIO_SIZE % (CHANNELS*SAMPLE_SIZE)) AUDIO_SIZE--;
    double dlit = AUDIO_SIZE/(SAMPLE_RATE*CHANNELS*SAMPLE_SIZE);
    printf("Duration: %.3lfs\n", dlit);

    return itype;
}

// =========================================================================================================================

int main(int argc, char *argv[])
{
    if (argc < 2) {
        puts("Invalid number of arguments\n");
        return 1;
    }

    std::string fname = argv[1];
    double mp = (argc > 2) ? atof(argv[2]) : 1;
    unsigned int FFT_SIZE = (argc > 3) ? atoi(argv[3]) : 2048;
    unsigned int SHIFT_SIZE = (mp == 1) ? FFT_SIZE : (FFT_SIZE>>1);           // if mp=1, it is recommended to shift full FFT_SIZE vector (not half) to avoid volume loss

    // DELTA_T_SIZE = (int)(DELTA_T*SAMPLE_RATE);          // AMPLITUDES
    // FFT_SIZE = 1 << ((int)(log2(DELTA_T_SIZE)) + 1);    // AMPLITUDES
    // DELTA_T_SIZE = 2205 amplitudes, if DELTA_T == 50 ms && SAMPLE_RATE == 44100 Hz
    // dj_fft library do not work correctly on audio data if vector size is 4096 or bigger
    // so manual set FFT_SIZE to 2048 (recommended)

    fp = fopen(fname.c_str(), "r");
    if (fp == NULL) {
        puts("File not found!\n");
        return 2;
    }

    if (!readHeader()) { puts("Cannot parse header.\n"); return 3; }

    // read data bytes
    unsigned char *adata = (unsigned char *) malloc(AUDIO_SIZE);
    fread(adata, 1, AUDIO_SIZE, fp);

    // save header
    char hdr[HDR_SIZE];
    fseek(fp, 0, 0);
    fread(hdr, 1, HDR_SIZE, fp);
    fclose(fp);

    //for (unsigned int i = 0; i != 64; ++i) printf("%u ", adata[i]); printf("\n");
    //for (unsigned int i = AUDIO_SIZE - 64; i != AUDIO_SIZE; ++i) printf("%u ", adata[i]); printf("\n");

    // process audio data
    std::vector<std::complex<float>> fftData[CHANNELS];
    unsigned int signedDiff = 1<<(8*SAMPLE_SIZE - 1);
    unsigned int i, l;
    unsigned short j, k;

    for (i = 0; i < AUDIO_SIZE; i += CHANNELS*SAMPLE_SIZE) {
        // bytes -> amplitudes
        for (j = 0; j != CHANNELS; ++j) {
            int smpl = 0;
            unsigned int i1 = i + j*SAMPLE_SIZE;
            if (BYTES_ORDER) memcpy(&smpl, adata+i1, SAMPLE_SIZE);      // big endian
            else {                                                      // little endian
                for (k = 0; k != SAMPLE_SIZE; ++k)      smpl |= (adata[i1 + k] << (8*k));
                if (smpl > signedDiff && !BYTE_FORMAT)  smpl -= (signedDiff << 1);
            }
            fftData[j].push_back(smpl);
        }

        // processing each chunk of FFT_SIZE size (except "tail" lesser than FFT_SIZE
        if (fftData[0].size() == FFT_SIZE) {
            for (j = 0; j != CHANNELS; ++j) {
                fftData[j] = fft1d(fftData[j], 1);      // fft forward (amplitudes -> frequencies)

                // pitch shift
                std::vector<std::complex<float>> fftCopy(FFT_SIZE);
                fftCopy[0] = fftData[j][0];
                for (l = 1; lround(l*mp) < SHIFT_SIZE; l++) fftCopy[lround(l*mp)] = fftData[j][l];  // (<=) mb?

                fftCopy = fft1d(fftCopy, -1);           // fft backward (frequencies -> amplitudes)
                fftData[j].clear();                     // reset fftData[]

                // amplitudes -> bytes (rewrite inidata)
                for (l = 0; l != FFT_SIZE; ++l)  {
                    int smpl = ROUND_2_INT(fftCopy[FFT_SIZE-1-l].real());
                    unsigned int bi = i - SAMPLE_SIZE*(l*CHANNELS - j);
                    for (k = 0; k != SAMPLE_SIZE; ++k) {
                        if (BYTES_ORDER)    adata[bi + SAMPLE_SIZE-1-k] = smpl & 0xff;    // big endian
                        else                adata[bi + k] = smpl & 0xff;                  // little endian
                        smpl >>= 8;
                    }
                }
            }
        }
    }

    // process tail
    i -= CHANNELS*SAMPLE_SIZE;
    if (fftData[0].size() > 0) {
        unsigned int lim = fftData[0].size();
        for (j = 0; j != CHANNELS; ++j) {
            while (fftData[j].size() != FFT_SIZE) fftData[j].push_back(0.f);

            fftData[j] = fft1d(fftData[j], 1);      // fft forward (amplitudes -> frequencies)

            // pitch shift
            std::vector<std::complex<float>> fftCopy(FFT_SIZE);
            fftCopy[0] = fftData[j][0];
            for (l = 1; lround(l*mp) < SHIFT_SIZE; l++) fftCopy[lround(l*mp)] = fftData[j][l];  // (<=) mb?

            fftCopy = fft1d(fftCopy, -1);           // fft backward (frequencies -> amplitudes)
            fftData[j].clear();                     // reset fftData[]

            // amplitudes -> bytes (rewrite inidata)
            for (l = 0; l != lim; ++l)  {
                int smpl = ROUND_2_INT(fftCopy[lim-1-l].real());
                unsigned int bi = i - SAMPLE_SIZE*(l*CHANNELS - j);
                for (k = 0; k != SAMPLE_SIZE; ++k) {
                    if (BYTES_ORDER)    adata[bi + SAMPLE_SIZE-1-k] = smpl & 0xff;    // big endian
                    else                adata[bi + k] = smpl & 0xff;                  // little endian
                    smpl >>= 8;
                }
            }
        }
    }

    //for (i = 0; i != 64; ++i) printf("%u ", adata[i]); printf("\n");
    //for (i = AUDIO_SIZE - 64; i != AUDIO_SIZE; ++i) printf("%u ", adata[i]); printf("\n");

    // save output
    fname.insert(fname.size()-4, "_" + std::to_string((int)(mp*100)) + "%");
    fp = fopen(fname.c_str(), "w");
    fwrite(hdr, 1, HDR_SIZE, fp);       // header
    fwrite(adata, 1, AUDIO_SIZE, fp);   // audio data

    free(adata);
    fclose(fp);

    return 0;
}

