%% RUN_unit_tests.m
% Lightweight unit tests for quasi-1D solver plumbing (GHOST-based BCs).
% Place in project root and run:
%   RUN_unit_tests

clc; clear; close all;

FUN_setup_path;

fprintf('\n====================================\n');
fprintf(' RUNNING UNIT TESTS (quasi-1D solver)\n');
fprintf('====================================\n');

try
    % ------------------------------------------------------------
    % Test 0: Initialize / finalize lifecycle
    % ------------------------------------------------------------
    fprintf('\n[Test 0] FUN_initialize()...\n');
    FUN_initialize();

    global MESH SOL BC GHOST
    assert(isstruct(MESH) && isfield(MESH,'N') && MESH.N >= 3, 'MESH not initialized correctly.');
    assert(isfield(MESH,'ng') && MESH.ng >= 1, 'MESH.ng missing or invalid.');

    assert(isstruct(SOL) && isfield(SOL,'U') && isfield(SOL,'Nv'), 'SOL not initialized correctly.');
    assert(isstruct(BC)  && isfield(BC,'left') && isfield(BC,'right'), 'BC not initialized correctly.');

    assert(isstruct(GHOST) && isfield(GHOST,'left') && isfield(GHOST,'right'), 'GHOST not initialized correctly.');
    assert(isfield(GHOST,'ng'), 'GHOST.ng missing.');

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
    % Test 2: Apply BC updates ghost states (no Uext)
    % ------------------------------------------------------------
    fprintf('\n[Test 2] FUN_apply_BC() updates GHOST.*.U sizes & finiteness...\n');

    % Snapshot before (to check something changes; not always guaranteed,
    % but useful during development)
    hasLeft  = isfield(GHOST.left,'U')  && ~isempty(GHOST.left.U);
    hasRight = isfield(GHOST.right,'U') && ~isempty(GHOST.right.U);
    if hasLeft,  Ugl_before = GHOST.left.U;  else, Ugl_before = []; end
    if hasRight, Ugr_before = GHOST.right.U; else, Ugr_before = []; end

    FUN_apply_BC();

    ng = MESH.ng;

    assert(isfield(GHOST.left,'U')  && ~isempty(GHOST.left.U),  'GHOST.left.U was not set by FUN_apply_BC().');
    assert(isfield(GHOST.right,'U') && ~isempty(GHOST.right.U), 'GHOST.right.U was not set by FUN_apply_BC().');

    assert(isequal(size(GHOST.left.U),  [ng, Nv]), 'GHOST.left.U must be ng x Nv.');
    assert(isequal(size(GHOST.right.U), [ng, Nv]), 'GHOST.right.U must be ng x Nv.');

    assert(all(isfinite(GHOST.left.U(:))),  'GHOST.left.U contains non-finite values.');
    assert(all(isfinite(GHOST.right.U(:))), 'GHOST.right.U contains non-finite values.');

    % Optional: check something actually changed (skip if BC is strictly transmissive)
    if ~isempty(Ugl_before)
        if max(abs(GHOST.left.U(:) - Ugl_before(:))) == 0
            fprintf('  (note) left ghost unchanged (may be OK for transmissive BC)\n');
        end
    end
    if ~isempty(Ugr_before)
        if max(abs(GHOST.right.U(:) - Ugr_before(:))) == 0
            fprintf('  (note) right ghost unchanged (may be OK for transmissive BC)\n');
        end
    end

    fprintf('  OK: GHOST.left/right.U are %d x %d and finite\n', ng, Nv);

    % ------------------------------------------------------------
    % Test 3: RHS call works + output dimensions/finiteness
    % ------------------------------------------------------------
    fprintf('\n[Test 3] FUN_get_RHS() basic sanity...\n');

    % Many codes compute RHS using SOL+GHOST globals, so no inputs.
    % Support both possible signatures:
    try
        out = cell(1,2);
        [out{:}] = FUN_get_RHS();  % try [A,RHS]
        A   = out{1};
        RHS = out{2};
    catch
        % try RHS-only
        A = [];
        RHS = FUN_get_RHS();
    end

    assert(isempty(A), 'A is expected to be empty for explicit stage.');
    assert(isequal(size(RHS), [Nv, N]), 'RHS must be Nv x N.');
    assert(all(isfinite(RHS(:))), 'RHS contains non-finite values.');

    fprintf('  OK: RHS is %d x %d\n', size(RHS,1), size(RHS,2));

    % ------------------------------------------------------------
    % Test 4 (optional): uniform state => RHS ~ 0
    % Works if BC + sources are uniform-consistent.
    % ------------------------------------------------------------
    fprintf('\n[Test 4] Uniform state => RHS near zero (optional)...\n');
    fprintf('\n Continue with this test? \nSelection of boundary values and source terms can lead to non zero RHS');
    fprintf('\n This does not neccessarily indicate wrong solver implementation');
    user_in = input('\n Continue? \t y/n','s');

    if strcmp(user_in,'y')
        U_save = SOL.U;

        SOL.U(:,:) = repmat(SOL.U(1,:), N, 1);
        FUN_apply_BC();  % must refresh ghosts for uniform state

        % RHS again
        try
            out = cell(1,2);
            [out{:}] = FUN_get_RHS();
            RHS_uni = out{2};
        catch
            RHS_uni = FUN_get_RHS();
        end

        rhs_norm_inf = max(abs(RHS_uni(:)));
        tol = 1e-10;  % relax during early dev if needed

        assert(rhs_norm_inf < tol, ...
            'Uniform state should give near-zero RHS. max|RHS|=%.3e', rhs_norm_inf);

        SOL.U = U_save;  % restore
        FUN_apply_BC();  % restore ghosts to match restored SOL

        fprintf('  OK: max|RHS| = %.3e < %.1e\n', rhs_norm_inf, tol);
    
    end

    fprintf('\n====================================\n');
    fprintf(' ALL UNIT TESTS PASSED \n');
    fprintf('====================================\n');

catch ME
    fprintf('\n====================================\n');
    fprintf(' UNIT TEST FAILED \n');
    fprintf('====================================\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);

    rethrow(ME);
end

% Always finalize
fprintf('\n[Cleanup] FUN_finalize()...\n');
try
    FUN_finalize();
catch
end
fprintf('Done.\n');
