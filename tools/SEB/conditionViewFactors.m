function [vf, svf, repaired] = conditionViewFactors(vf, areas)
% Validate or conservatively repair open-domain View3D factors.
% This mirrors the default Python preprocessing policy.
tolerance = 1e-3;
closure_margin = 5e-4;
max_row_sum = 2.0;
max_overfull_area_fraction = 0.05;
max_exchange_reduction_fraction = 0.01;
max_reciprocity_error = 0.01;

vf = sparse(double(vf));
areas = double(areas(:));
nfacets = size(vf, 1);
if size(vf, 2) ~= nfacets || numel(areas) ~= nfacets
    error('View-factor matrix and facet-area dimensions do not match.')
end
if any(~isfinite(areas)) || any(areas <= 0)
    error('Facet areas must be finite and strictly positive.')
end

[rows, cols, values] = find(vf);
if any(~isfinite(values))
    error('View3D factors contain non-finite values.')
end
if any(values < -tolerance)
    error('View3D factors contain values below the physical tolerance.')
end
if any(values > 1 + tolerance)
    error('View3D factors contain individual values above the physical tolerance.')
end
tiny_negative = values < 0;
repaired = any(tiny_negative);
if repaired
    values(tiny_negative) = 0;
    vf = sparse(rows, cols, values, nfacets, nfacets);
end

exchange = spdiags(areas, 0, nfacets, nfacets) * vf;
reciprocity_before = reciprocity_error(exchange);
if reciprocity_before > max_reciprocity_error
    error('View3D reciprocity L1 error %.6g exceeds the automatic-repair limit %.6g.', ...
          reciprocity_before, max_reciprocity_error)
end

row_sum_raw = full(sum(vf, 2));
if max(row_sum_raw) > max_row_sum
    error('Maximum View3D row sum %.6g exceeds the automatic-repair limit %.6g.', ...
          max(row_sum_raw), max_row_sum)
end
if ~any(row_sum_raw > 1 + tolerance)
    svf = max(1 - row_sum_raw, 0);
    validate_conditioned_view_factors(vf, svf, areas, tolerance, max_reciprocity_error);
    return
end

inverse_counts = spones(exchange) + spones(exchange');
[count_rows, count_cols, count_values] = find(inverse_counts);
inverse_counts = sparse(count_rows, count_cols, 1 ./ count_values, nfacets, nfacets);
symmetric_exchange = (exchange + exchange') .* inverse_counts;
row_sum = full(sum(symmetric_exchange, 2)) ./ areas;
overfull = row_sum > 1;
materially_overfull = row_sum > 1 + tolerance;
if max(row_sum) > max_row_sum
    error('Reciprocity adjustment gives row sum %.6g above the automatic-repair limit %.6g.', ...
          max(row_sum), max_row_sum)
end

overfull_area_fraction = sum(areas(materially_overfull)) / sum(areas);
if overfull_area_fraction > max_overfull_area_fraction
    error('Materially overfull facets represent %.3f%% of area; limit is %.3f%%.', ...
          100 * overfull_area_fraction, 100 * max_overfull_area_fraction)
end

scale = ones(nfacets, 1);
scale(overfull) = (1 - closure_margin) ./ row_sum(overfull);
[rows, cols, values] = find(symmetric_exchange);
values = values .* min(scale(rows), scale(cols));
repaired_exchange = sparse(rows, cols, values, nfacets, nfacets);
exchange_before = sum(nonzeros(symmetric_exchange));
exchange_reduction_fraction = max(exchange_before - sum(nonzeros(repaired_exchange)), 0) / exchange_before;
if exchange_reduction_fraction > max_exchange_reduction_fraction
    error('Repair would remove %.3f%% of exchange area; limit is %.3f%%.', ...
          100 * exchange_reduction_fraction, 100 * max_exchange_reduction_fraction)
end

vf = spdiags(1 ./ areas, 0, nfacets, nfacets) * repaired_exchange;
row_sum_after = full(sum(vf, 2));
svf = max(1 - row_sum_after, 0);
validate_conditioned_view_factors(vf, svf, areas, tolerance, max_reciprocity_error);
repaired = true;
warning('uDALES:ViewFactorRepair', ...
        ['View3D open-domain closure repaired: %d materially overfull rows, ' ...
         'maximum row sum %.6g -> %.6g, materially overfull facet area %.3f%%, ' ...
         'removed exchange area %.3f%%.'], ...
        nnz(materially_overfull), max(row_sum), max(row_sum_after), ...
        100 * overfull_area_fraction, 100 * exchange_reduction_fraction)
end

function value = reciprocity_error(exchange)
denominator = sum(abs(nonzeros(exchange)));
if denominator == 0
    value = 0;
else
    value = sum(abs(nonzeros(exchange - exchange'))) / denominator;
end
end

function validate_conditioned_view_factors(vf, svf, areas, tolerance, reciprocity_tolerance)
closure = full(sum(vf, 2)) + svf;
if any(~isfinite(svf)) || any(svf < -tolerance) || any(svf > 1 + tolerance)
    error('Conditioned sky view factors are outside physical bounds.')
end
if any(abs(closure - 1) > tolerance)
    error('Conditioned view factors fail surface-plus-sky closure.')
end
if reciprocity_error(spdiags(areas, 0, numel(areas), numel(areas)) * vf) > reciprocity_tolerance
    error('Conditioned view factors fail reciprocity validation.')
end
end
