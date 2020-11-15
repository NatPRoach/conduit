import math
import sequtils


proc findLocalMinima(region : seq[float]) : seq[uint32] =
  for i in 0..<region.len:
    var previous, next = 0.0
    if i > 0:
      previous = region[i-1]
    if i < region.len - 1:
      next = region[i+1]
    if region[i] > previous and region[i] > next:
      result.add(uint32(i))


proc calculateKDEmaxima(weighted_positions : seq[(uint32,uint32)], calculation_window : uint16 = 30'u16, stddev : float = 10.0 ) : seq[uint32] =
  let start_pos = weighted_positions[0][0]
  let end_pos = weighted_positions[^1][0]
  let delta = end_pos - start_pos
  var region = repeat(0.0, delta + 1)
  for j, (pos, weight) in weighted_positions:
    let offset = int(pos) - int(start_pos)
    var low_region_idx = offset - int(calculation_window)
    var low_early_cutoff = 0
    var high_region_idx = offset + int(calculation_window)
    if low_region_idx < 0:
      low_early_cutoff = abs(low_region_idx)
      low_region_idx = 0
    if high_region_idx >= region.len:
      high_region_idx = region.len - 1

    for i in low_region_idx..high_region_idx:
      let numer = float(weight) * math.exp(-(math.pow((float(pos) - float(i)), 2)/(2 * math.pow(stddev, 2))))
      let denom = stddev * math.sqrt(2 * math.PI)
      region[i] += numer / denom
  let offset_maxima = findLocalMinima(region)
  for maxima in offset_maxima:
    result.add(start_pos + maxima)

proc getKDEmaxima*(weighted_positions : seq[(uint32,uint32)], grouping_window : uint16 = 15'u16, calculation_window : uint16 = 30'u16, stddev : float = 10.0 ) : seq[uint32] = 
  # For an array of genome positions and weights,
  # find local maxima of potential splice junctions based on weighted KDE
  var start_idx, end_idx = 0
  var last_position = weighted_positions[0][0]
  for i in 1..<weighted_positions.len:
    if weighted_positions[i][0] - last_position > grouping_window:
      for maxima in calculateKDEmaxima(weighted_positions[start_idx..end_idx],calculation_window,stddev):
        result.add(maxima)
      start_idx = i
    end_idx = i
    last_position = weighted_positions[i][0]
  for maxima in calculateKDEmaxima(weighted_positions[start_idx..end_idx],calculation_window,stddev):
    result.add(maxima)