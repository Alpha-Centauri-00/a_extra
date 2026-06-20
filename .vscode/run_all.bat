@echo off
setlocal

set "SCRIPT_DIR=src"

echo ==============================
echo Starting Galaxy Pipeline
echo ==============================

for %%S in (
    a_extract_morphology.py
    b_compute_J_visible.py
    c_fit_k_continuous.py
    d_merge_k_morphology_correlations.py
    e_lt_analysis_galaxies.py
    f_cluster_spin_coupling.py
    g_mond_comparison.py
    z_Galaxy_explorer_ml_sensitivity.py
) do (
    echo.
    echo Running %%S...
    uv run python "%SCRIPT_DIR%\%%S"
    if errorlevel 1 (
        echo ERROR: %%S failed!
        exit /b 1
    )
    echo Completed %%S
)

echo.
echo ==============================
echo Pipeline completed successfully
echo ==============================
exit /b 0