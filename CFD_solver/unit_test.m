%% RUN_unit_tests.m
% Lightweight unit tests for quasi-1D solver plumbing.
% Place in project root and run:
%   RUN_unit_tests

clc; clear; close all;

fprintf('\n====================================\n');
fprintf(' RUNNING UNIT TESTS (quasi-1D solver)\n');
fprintf('====================================\n');

try
    % ------------------------------------------------------------
    % Test 0: Initialize / finalize lifecycle
    % ------------------------------------------------------------
    fprintf('\n[Test 0] FUN_initialize()...\n');
    FUN_initialize();

    global MESH SOL BC
    assert(isstruct(MESH) && isfield(MESH,'N') && MESH.N >= 3, 'MESH not initialized correctly.');
    assert(isfield(MESH,'ng') && isfield(MESH,'Ntot'), 'Ghost indexing fields missing in MESH.');
    assert(isfield(MESH,'iC') && isfield(MESH,'phys2ext'), 'MESH.iC or MESH.phys2ext missing.');

    assert(isstruct(SOL) && isfield(SOL,'U') && isfield(SOL,'Nv'), 'SOL not initialized correctly.');
    assert(isstruct(BC)  && isfield(BC,'left') && isfield(BC,'right'), 'BC not initialized correctly.');

    fprintf('  OK\n');

    % ------------------------------------------------------------
    % Test 1: Solution dimensions and finiteness
    % ------------------------------------------------------------
    fprintf('\n[Test 1] SOL.U dimensions / finiteness...\n');
    N  = MESH.N;
    Nv = SOL.Nv;

    assert(isequal(size(SOL.U), [N, Nv]), 'SOL.U must be N x Nv.');
    assert(all(isfinite(SOL.U(:))), 'SOL.U contains non-finite values.');

    fprintf('  OK: SOL.U is %d x %d\n', N, Nv);

    % ------------------------------------------------------------
    % Test 2: Apply BC -> Uext dimensions and embedding
    % ------------------------------------------------------------
    fprintf('\n[Test 2] FUN_apply_BC() -> Uext size & embedding...\n');

    % If your BC function has a different name, change it here:
    Uext = FUN_apply_BC();

    ng = MESH.ng;
    assert(isequal(size(Uext), [N+2*ng, Nv]), 'Uext must be (N+2ng) x Nv.');

    % Physical core in Uext must match SOL.U exactly
    Ucore = Uext(MESH.iC, :);
    assert(max(abs(Ucore(:) - SOL.U(:))) == 0, 'Uext(MESH.iC,:) must equal SOL.U.');

    fprintf('  OK: Uext is %d x %d, core matches SOL.U\n', size(Uext,1), size(Uext,2));

    % ------------------------------------------------------------
    % Test 3: RHS call works + output dimensions/finiteness
    % ------------------------------------------------------------
    fprintf('\n[Test 3] FUN_get_RHS(Uext) basic sanity...\n');
    [A, RHS] = FUN_get_RHS(Uext); %#ok<ASGLU>  % A is empty for now

    assert(isempty(A), 'A is expected to be empty for explicit stage.');
    assert(isequal(size(RHS), [N, Nv]), 'RHS must be N x Nv.');
    assert(all(isfinite(RHS(:))), 'RHS contains non-finite values.');

    fprintf('  OK: RHS is %d x %d\n', size(RHS,1), size(RHS,2));

    % ------------------------------------------------------------
    % Test 4 (optional but very useful): uniform state => RHS ~ 0
    % This assumes transmissive BC and zero (or uniform-consistent) sources.
    % ------------------------------------------------------------
    fprintf('\n[Test 4] Uniform state => RHS near zero (optional)...\n');

    % Save current state
    U_save = SOL.U;

    % Force uniform state based on first cell
    SOL.U(:,:) = repmat(SOL.U(1,:), N, 1);

    Uext = FUN_apply_BC();
    [~, RHS_uni] = FUN_get_RHS(Uext);

    rhs_norm_inf = max(abs(RHS_uni(:)));
    tol = 1e-10;  % relax if needed during early development
    assert(rhs_norm_inf < tol, 'Uniform state should give near-zero RHS. max|RHS|=%.3e', rhs_norm_inf);

    % Restore
    SOL.U = U_save;

    fprintf('  OK: max|RHS| = %.3e < %.1e\n', rhs_norm_inf, tol);

    fprintf('\n====================================\n');
    fprintf(' ALL UNIT TESTS PASSED \n');
    fprintf('====================================\n');

catch ME
    fprintf('\n====================================\n');
    fprintf(' UNIT TEST FAILED \n');
    fprintf('====================================\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);

    % Re-throw after cleanup so you can click into the error in MATLAB
    rethrow(ME);
end

% Always finalize (even if a test failed)
fprintf('\n[Cleanup] FUN_finalize()...\n');
try
    FUN_finalize();
catch
    % ignore finalize errors
end
fprintf('Done.\n');
