# Run profiles

A profile is a named set of settings. `PROFILE=<name>` loads `config/profiles/<name>.env`.

Precedence, highest first:

    environment  >  profile  >  registry default

so a profile fixes the shape of a run while a single environment variable still overrides one
setting for a one-off. The run log and `run_manifest.json` record which layer supplied each value,
so "what did this run use" never depends on remembering what was set where.

Profiles are version-controlled on purpose: two runs are compared by diffing their profiles, and a
profile named in a manifest is recoverable months later.

Keys are the `name` column of `config/knobs.json`. A key not in the registry is rejected rather than
ignored, because a silently-ignored setting is the failure this whole mechanism exists to stop.

**Model vs operational.** Keys marked `model` in the registry change what is estimated; the defaults
are the validated standard (see `docs/sampler_model_specification.md`). A profile that sets one is
making a modelling claim, and the run announces it as a SPEC DEVIATION. Operational keys — cores,
sweeps, paths — never affect results.
