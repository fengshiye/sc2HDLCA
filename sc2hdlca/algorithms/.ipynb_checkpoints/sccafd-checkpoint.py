import subprocess
import os
import importlib.resources

import subprocess
from pathlib import Path
import importlib.resources

def run_SCCAF_D(bulk=None, python_home=None, results_path=None,R_home=None):
    # 找到包里的 R 脚本目录
    r_dir = importlib.resources.files("sc2hdlca").joinpath("R/SCCAF-D")
    r_script = r_dir / "run_SCCAF_D.R"

    # 用 cwd 确保 R 的 source 识别
    subprocess.run(
        [R_home, str(r_script), bulk, python_home, results_path],
        check=True,
        cwd=str(r_dir)
    )