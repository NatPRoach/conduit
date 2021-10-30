import math
import sequtils


proc findLocalMinima(region: seq[float]): seq[uint32] =
  ##
  for i in 0..<region.len:
    var previous, next = 0.0
    if i > 0:
      previous = region[i-1]
    if i < region.len - 1:
      next = region[i+1]
    if region[i] > previous and region[i] > next:
      result.add(uint32(i))


proc calculateKDEmaxima(weighted_positions: seq[(uint32, uint32)],
                        calculation_window: uint16 = 30'u16,
                        stddev: float = 10.0): seq[uint32] =
  ##
  let startPos = weighted_positions[0][0]
  let endPos = weighted_positions[^1][0]
  let delta = endPos - startPos
  var region = repeat(0.0, delta + 1)
  for j, (pos, weight) in weighted_positions:
    let offset = int(pos) - int(startPos)
    var lowRegionIdx = offset - int(calculation_window)
    var lowEarlyCutoff = 0
    var highRegionIdx = offset + int(calculation_window)
    if lowRegionIdx < 0:
      lowEarlyCutoff = abs(lowRegionIdx)
      lowRegionIdx = 0
    if highRegionIdx >= region.len:
      highRegionIdx = region.len - 1

    for i in lowRegionIdx..highRegionIdx:
      let numer = float(weight) *
                  math.exp(-(math.pow((float(pos) - float(i)), 2) /
                  (2 * math.pow(stddev, 2))))
      let denom = stddev * math.sqrt(2 * math.PI)
      region[i] += numer / denom
  let offsetMaxima = findLocalMinima(region)
  for maxima in offsetMaxima:
    result.add(startPos + maxima)

proc getKDEmaxima*(weighted_positions: seq[(uint32, uint32)],
                   grouping_window: uint16 = 15'u16,
                   calculation_window: uint16 = 30'u16,
                   stddev: float = 10.0): seq[uint32] =
  ## For an array of genome positions and weights,
  ## find local maxima of potential splice junctions based on weighted KDE
  var startIdx, endIdx = 0
  var lastPosition = weighted_positions[0][0]
  for i in 1..<weighted_positions.len:
    if weighted_positions[i][0] - lastPosition > grouping_window:
      for maxima in calculateKDEmaxima(weighted_positions[startIdx..endIdx],
                                       calculation_window,
                                       stddev):
        result.add(maxima)
      startIdx = i
    endIdx = i
    lastPosition = weighted_positions[i][0]
  for maxima in calculateKDEmaxima(weighted_positions[startIdx..endIdx],
                                   calculation_window,
                                   stddev):
    result.add(maxima)
