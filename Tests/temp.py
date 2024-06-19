def sl(start, end, step):
    result = []
    current = start
    while current != end:
        result.append(current)
        current += step
    return result