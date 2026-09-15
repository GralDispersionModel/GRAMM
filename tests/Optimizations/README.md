# Optimization regression checks

Build the original V2701 core in a separate checkout and build this branch into a different directory. Use .NET 10. Run each harness against the original and patched assemblies, then compare the result JSON hashes.

```powershell
dotnet build src/GRAMM.csproj -c Release -o artifacts/patched
dotnet build tests/Optimizations -c Release -o artifacts/harness
dotnet artifacts/harness/Optimizations.dll ../baseline/src/bin/Release/net10.0/GRAMM.dll artifacts/baseline.json
dotnet artifacts/harness/Optimizations.dll artifacts/patched/GRAMM.dll artifacts/patched.json
```

The harness loads production methods from the supplied assembly. Use a fresh process for each measurement. Timings describe the changed stage and vary by machine; they are not total-model speedup estimates.

The harness calls `Define_Arrays`, measures retained managed memory, initializes a reproducible pressure-correction fixture and hashes all modified numerical fields. It also exercises array clearing.

For full runs, obtain the official GUI repository and use its `SampleProjects/AscendingBridge/Computation` directory as read-only input. Run in a real Windows console; the existing core uses console cursor APIs and cannot run this case with stdout redirected to a pipe.

```powershell
py -3 tests/Optimizations/validation.py --baseline ../baseline/src/bin/Release/net10.0/GRAMM.dll --patched artifacts/patched/GRAMM.dll --sample-computation ../GUI/SampleProjects/AscendingBridge/Computation --output artifacts/full-regression
```

The script copies an input allowlist into a new output directory, uses one worker and short integration times, and compares default plus custom-temperature situations. It compares the contents of scalar ZIP containers, so archive timestamps do not count as numerical changes. Python uses only its standard library.
