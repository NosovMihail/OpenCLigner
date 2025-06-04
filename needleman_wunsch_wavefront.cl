__kernel void needleman_wunsch_wavefront(
    __global const char* seq1,
    __global const char* seq2,
    __global int* matrix,
    __global int* max_score,
    const int len1,
    const int len2,
    const int diag,
    const int match_score,
    const int mismatch_penalty,
    const int gap_penalty,
    __global int* max_i,
    __global int* max_j)
{
    int gid = get_global_id(0);
    int i = gid + 1;
    int j = diag - i + 1;

    if (i <= 0 || j <= 0 || i > len1 || j > len2)
        return;

    char a = seq1[i - 1];
    char b = seq2[j - 1];
    int match = (a == b) ? match_score : mismatch_penalty;

    int idx     = i * (len2 + 1) + j;
    int up      = matrix[(i - 1) * (len2 + 1) + j] + gap_penalty;
    int left    = matrix[i * (len2 + 1) + (j - 1)] + gap_penalty;
    int diagval = matrix[(i - 1) * (len2 + 1) + (j - 1)] + match;

    int score = max(diagval, max(up, left));
    matrix[idx] = score;
    *max_i = len1;
    *max_j = len2;
}
