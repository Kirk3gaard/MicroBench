{
  description = "A flake for the trustmap repo";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        rEnv = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            tidyverse
            languageserver
            janitor
            data_table
            readxl
            kableExtra
          ];
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            rEnv
            # For file watching and plotting workflow
            pkgs.inotify-tools
            pkgs.feh
            pkgs.xdg-utils
            # For building
            pkgs.openssl
            pkgs.pkg-config

            pkgs.libpng
          ];

          shellHook = ''
            export R_LIBS_USER=""
          '';
        };
      }
    );
}
