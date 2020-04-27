function PSNR_val = psnr(ref, x)
max_val = max(ref(:));
ref = ref / max_val;
x = x / max_val;
PSNR_val = 10 * log10(1/mean((ref(:) - x(:)).^2));
