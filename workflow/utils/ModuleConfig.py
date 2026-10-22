from pathlib import Path
from glob import glob
import os
import yaml
import copy

def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge two dictionaries (override wins)."""
    for k, v in override.items():
        if isinstance(v, dict) and k in base:
            deep_merge(base[k], v)
        else:
            base[k] = v
    return base

class ModuleConfig:
    """
    Handles module parameters, paths, and input detection.
    """

    def __init__(self, module_name: str, global_config: dict):
        self.module_name = module_name
        self.global_config = global_config

        # Load default config
        self.config = global_config
        with open(Path(__file__).parents[1] / module_name / "config.defaults.yaml") as f:
            defaults = yaml.safe_load(f) or {}
        module_defaults = defaults.get(module_name, {})
        user_cfg = global_config.get(module_name, {})
        self.config = deep_merge(copy.deepcopy(module_defaults), user_cfg)

        # Output directories
        output_dir = Path(global_config.get("output_dir", "results"))
        self.temp_dir = output_dir / "data" / "temp" / self.module_name
        self.out_dir = output_dir / "data" / self.module_name
        self.reports_dir = output_dir / "reports"
        self.bench_dir = output_dir / "benchmarking" / self.module_name
        self.log_dir = output_dir / "logs" / self.module_name

        # Define Inputs, Conda envs
        self.set_inputs()
        self.set_env_paths()

    def get(self, key, default=None):
        return self.config.get(key, default)
    
    def set_inputs(self):
        """
        Determine input files for the module based on the configuration.

        Handles three cases:
        1. Input is specified as another module name.
        2. Input is a single file path ending with '.h5ad' or '.zarr'.
        3. Input is a directory containing '.h5ad' or '.zarr' files.

        Sets:
        - self.input_dir: directory where input files are located.
        - self.input_files: list of input file names without extensions.
        - self.input_type: type of input, either 'h5ad' or 'zarr', used for Snakemake rules.
        """
        input_map = self.global_config.get("input", {})
        input_str = input_map.get(self.module_name, None)
        self.input_files = []
        self.input_type = None

        # Case 1: Input string refers to another module
        if input_str in self.global_config["modules"]:
            self.input_type = "zarr"
            self.input_dir = Path(self.global_config["output_dir"]) / "data" / input_str
            
            input_config = ModuleConfig(input_str, self.global_config)
            self.input_files = [f for f in input_config.input_files]
            
            # add preprocessing combos to input files names if using preprocessing outputs
            if input_str == "preprocessing":
                self.input_files = [
                    f"{f}.{combo}"
                    for f in self.input_files
                    for combo in self.config["prepro_combos"]
                ]

        # Case 2: Input is a single file OR a single zarr store
        elif os.path.isfile(input_str) or (
            os.path.isdir(input_str) and input_str.endswith(".zarr")
        ):
            if input_str.endswith('.h5ad'):
                self.input_type = 'h5ad'
            elif input_str.endswith('.zarr'):
                self.input_type = 'zarr'
            else:
                raise ValueError(f"Unsupported file type for input: {input_str}")

            self.input_dir = os.path.dirname(input_str)
            self.input_files = [os.path.splitext(os.path.basename(input_str))[0]]

        # Case 3: Input is a directory of files
        elif os.path.isdir(input_str):
            self.input_dir = input_str
            for f in os.listdir(input_str):
                if f.endswith('.h5ad'):
                    self.input_type = 'h5ad'
                    self.input_files.append(os.path.splitext(f)[0])
                elif f.endswith('.zarr'):
                    self.input_type = 'zarr'
                    self.input_files.append(os.path.splitext(f)[0])

    def set_env_paths(self):
        """
        Set up environment YAML paths for the module and related tools.

        Returns a dictionary `self.env_paths` with:
        - Module-specific environment YAML.
        - Plotting environment YAML.
        - Any tool-specific environments defined in the module's configuration.

        The paths are relative to the project root (two levels up from this file).
        """
        PROJECT_ROOT = Path(__file__).resolve().parents[2]
        self.env_paths = {}

        # Set module-wide env path
        self.env_paths[self.module_name] = PROJECT_ROOT / "envs" / f"{self.module_name}.yaml"

        # Set plotting specific env paths
        self.env_paths['reports'] = PROJECT_ROOT / "envs" / "reports.yaml"

        # Add a tool.yaml env if tool defined in config
        for submodule_name, submodule_cfg in self.config.items():
            if isinstance(submodule_cfg, dict):
                if "tool" in submodule_cfg:
                    tools = submodule_cfg["tool"]
                    if isinstance(tools, str):
                        tools = [tools]
                    for tool in tools:
                        self.env_paths[tool] = PROJECT_ROOT / "envs" / f"{tool}.yaml"
        
        # Add path to integration envs dir
        if self.module_name == "preprocessing":
            self.env_paths['integration'] = PROJECT_ROOT / "envs" / "integration"

    # QC MODULE SPECIFIC FUNCTIONS
    def get_doublet_tools(self):
        """Return list of doublet detection tools."""
        tools = self.config.get("doublets", {}).get("tool", [])
        if not isinstance(tools, list):
            tools = [tools]
        return tools
    
    # PREPROCESSING MODULE SPECIFIC FUNCTIONS
    def get_gene_lists_path(self, cc_method):
        """
        Retrieve full paths to gene list files used for cell cycle analysis.

        Looks for:
        - S-phase genes file
        - G2M-phase genes file

        Returns:
            List[Path]: List of full paths to the gene list files or empty list.
        """
        files = []

        if cc_method == "scanpy":
            PROJECT_ROOT = Path(__file__).resolve().parents[2]
            for l in ["s_genes_file","g2m_genes_file"]:
                filename = self.config["cell_phase_scoring"][l]
                filepath = PROJECT_ROOT / 'resources' / filename
                files.append(filepath)
        
        return files

    def get_layer_key(self, integration_method):
        if integration_method in ["mnn","combat"]:
            return f"{integration_method}_counts"
        elif integration_method in ["shuffle-per-batch", "shuffle-per-label"]:
            return f"{integration_method}_counts"
        else:
            return "norm_counts"
    
    def get_embedding_key(self, integration_method):
        if integration_method in [
            "harmony","scanorama","css","cca","rpca","scvi","scanvi", "concord", #PCA embedding modifying/replacing methods
            "perfect-embed-label", "perfect-embed-label-jitter", "batchwise-pca", #PCA embedding modifying/replacing control methods 
            ]: 
            embed_key = f"X_{integration_method}"
        else:
            embed_key = "X_pca"
        return embed_key
    
    def get_integration_product(self, integration_method):
        """
        Returns the integrated/corrected space produced by integration method.
        """
        step_produces = {
            # classic steps
            "pca":       "embedding_space",
            "neighbors": "graph_space",
            "projection":None,
            "noint":None,
            # integration methods
            "combat":    "feature_space",
            "mnn":       "feature_space",
            "scvi":      "embedding_space",
            "scanvi":    "embedding_space",
            "harmony":   "embedding_space",
            "scanorama": "embedding_space",
            "css":       "embedding_space",
            "cca":       "embedding_space",
            "rpca":      "embedding_space",
            "bbknn":     "graph_space",
            "concord":   "embedding_space",
            # integration controls
            "shuffle":   "feature_space",
            "shuffle-per-batch":   "feature_space",
            "shuffle-per-label":   "feature_space",
            "batchwise-pca":      "embedding_space",
            "perfect-embed-label":        "embedding_space",
            "perfect-embed-label-jitter": "embedding_space",
        }

        return step_produces[integration_method]

    def get_prepro_input(self, step_name, wildcards):
        """
        Return the single required input file for a preprocessing/integration step.
        """

        step_order = {
            # classic steps
            "noint":    ["pca", "neighbors", "projection"],
            # integration steps
            "combat":   ["combat", "pca", "neighbors", "projection"],
            "mnn":      ["mnn", "pca", "neighbors", "projection"],
            "harmony":  ["pca", "harmony", "neighbors", "projection"],
            "scanorama":["pca", "scanorama", "neighbors", "projection"],
            "scvi":     ["scvi", "neighbors", "projection"],
            "scanvi":   ["scanvi", "neighbors", "projection"],
            "bbknn":    ["pca", "bbknn", "projection"],
            "concord":  ["concord", "neighbors", "projection"],
            "css":      ["pca", "css", "neighbors", "projection"],
            "cca":      ["pca", "cca", "neighbors", "projection"],
            "rpca":     ["pca", "rpca", "neighbors", "projection"],
            # integration controls
            "shuffle":      ["shuffle", "pca", "neighbors", "projection"],
            "shuffle-per-batch":      ["shuffle-per-batch", "pca", "neighbors", "projection"],
            "shuffle-per-label":      ["shuffle-per-label", "pca", "neighbors", "projection"],
            "batchwise-pca":         ["batchwise-pca", "neighbors", "projection"],
            "perfect-embed-label":           ["perfect-embed-label", "neighbors", "projection"],
            "perfect-embed-label-jitter":    ["perfect-embed-label-jitter", "neighbors", "projection"],
        }

        step_consumes = {
            # classic steps
            "pca":       "feature_space",
            "neighbors": "embedding_space",
            "projection":"graph_space",
            # integration steps
            "combat":    "feature_space",
            "mnn":       "feature_space",
            "scvi":      "feature_space",
            "scanvi":    "feature_space",
            "concord":   "feature_space",
            "harmony":   "embedding_space",
            "scanorama": "embedding_space",
            "css":       "embedding_space",
            "cca":       "embedding_space",
            "rpca":      "embedding_space",
            "bbknn":     "embedding_space",
            # integration control steps
            "shuffle":              "feature_space",
            "shuffle-per-batch":     "feature_space",
            "shuffle-per-label":     "feature_space",
            "batchwise-pca":        "feature_space",
            "perfect-embed-label":          "feature_space",
            "perfect-embed-label-jitter":   "feature_space",
        }

        folder_seed = Path(self.temp_dir) / wildcards.filename
        integration = wildcards.integration
        feature_space = wildcards.feature_space
        chain = step_order[integration]

        # Initial state (no integration)
        state = {
            "feature_space": (Path(f"{folder_seed}.feature_selection.zarr")),
            "embedding_space": (Path(f"{folder_seed}.{feature_space}_noint.embedding_space.pca.zarr")),
            "graph_space": (Path(f"{folder_seed}.{feature_space}_noint.graph_space.neighbors.zarr")),
        }

        # Walk the chain and update state for integration methods
        if wildcards.integration != "noint":
            for step in chain:
                if step == step_name: #until we reach step_name
                    break
                if chain.index(step_name) > chain.index(wildcards.integration): #do not update steps before integration
                    if self.get_integration_product(step) == "feature_space":
                        state["feature_space"] = (Path(f"{folder_seed}.{feature_space}_{integration}.feature_space.{step}.zarr"))
                    elif self.get_integration_product(step) == "embedding_space":
                        state["embedding_space"] = (Path(f"{folder_seed}.{feature_space}_{integration}.embedding_space.{step}.zarr"))
                    elif self.get_integration_product(step) == "graph_space":
                        state["graph_space"] = (Path(f"{folder_seed}.{feature_space}_{integration}.graph_space.{step}.zarr"))

        consumes = step_consumes[step_name]

        return state[consumes]

    def get_post_r_input(self, space, wildcards):
        # integration methods outputs are saved in npz for easy loading in python
        folder_seed = Path(self.temp_dir) / wildcards.filename
        integration = wildcards.integration
        feature_space = wildcards.feature_space
        state = {
            "feature_space": (Path(f"{folder_seed}.{feature_space}_{integration}.feature_selection.{integration}.npz")),
            "embedding_space": (Path(f"{folder_seed}.{feature_space}_{integration}.embedding_space.{integration}.npz")),
            "graph_space": (Path(f"{folder_seed}.{feature_space}_{integration}.graph_space.{integration}.npz")),
        }
        return state[space]
    