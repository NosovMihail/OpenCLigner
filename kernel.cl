__kernel void smith_waterman(
    __global const char* seq1,
    __global const char* seq2,
    __global int* matrix,
    __global int* max_score,
    const int len1,
    const int len2,
    const int match_score,
    const int mismatch_penalty,
    const int gap_penalty)
{
    int i = get_global_id(0) + 1;
    int j = get_global_id(1) + 1;

    if (i > len1 || j > len2) return;

    char a = seq1[i - 1];
    char b = seq2[j - 1];
    int match = (a == b) ? match_score : mismatch_penalty;

    int idx = i * (len2 + 1) + j;
    int up   = matrix[(i - 1) * (len2 + 1) + j] + gap_penalty;
    int left = matrix[i * (len2 + 1) + (j - 1)] + gap_penalty;
    int diag = matrix[(i - 1) * (len2 + 1) + (j - 1)] + match;

    int score = max(0, max(diag, max(up, left)));
    matrix[idx] = score;

    atomic_max(max_score, score);
}
