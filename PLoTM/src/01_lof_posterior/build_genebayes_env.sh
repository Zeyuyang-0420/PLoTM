#!/bin/bash
# Build a dedicated GPU conda env for the paper's GeneBayes LoF model
# (torch + ngboost + torchquad + xgboost). CUDA 12.4 driver -> cu121 wheels.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_PREFIX="/mnt/scratch/ZY2/.envs/env_genebayes"

if [ ! -d "$ENV_PREFIX" ]; then
  conda create -y -p "$ENV_PREFIX" python=3.10
fi
conda activate "$ENV_PREFIX"

python -m pip install --upgrade pip
# torch built for CUDA 12.1 (runs on the 550/CUDA-12.4 driver)
python -m pip install torch --index-url https://download.pytorch.org/whl/cu121
# GeneBayes model deps
python -m pip install "numpy<2" pandas scipy scikit-learn xgboost ngboost torchquad openpyxl

echo "=== verify ==="
python - <<'PY'
import torch, ngboost, torchquad, xgboost, scipy, pandas, numpy
print("torch", torch.__version__, "cuda?", torch.cuda.is_available(),
      torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
print("ngboost", ngboost.__version__, "xgboost", xgboost.__version__,
      "numpy", numpy.__version__)
import torchquad; print("torchquad ok")
PY
echo "ENV_READY"
