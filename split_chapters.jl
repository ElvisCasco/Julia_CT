#!/usr/bin/env julia
# split_chapters.jl
# Read Cameron & Trivedi_stata_vJulia.qmd, generate per-chapter qmd files
# in the same folder. Each chapter file contains:
#   1. YAML header (julia engine)
#   2. The shared preamble (imports + compat shims) — always included
#   3. ONLY the chunk-1 helpers that the chapter (transitively) uses
#   4. The chapter's example chunks
#
# Helper extraction is heuristic: top-level `function NAME(...)`,
# `const NAME = ...`, and `struct NAME` blocks (no leading indent),
# with any preceding `"""..."""` docstring attached.

const QMD = "Cameron & Trivedi_stata_vJulia.qmd"
src_lines = readlines(QMD)
n_lines   = length(src_lines)

# ----------------------------------------------------------------------
# 1. Find chapter boundaries — lines like `# 1 Stata basics`.
#    The regex limits N to 1-2 digits to avoid `# 1000 draws...` etc.
# ----------------------------------------------------------------------
const CHAPTER_RE = r"^# (\d{1,2}) (.+)$"
chapters = NamedTuple{(:num, :title, :start_line),Tuple{Int,String,Int}}[]
for (i, line) in enumerate(src_lines)
    m = match(CHAPTER_RE, line)
    if m !== nothing
        n = parse(Int, m.captures[1])
        if 1 <= n <= 18 && (isempty(chapters) || n == chapters[end].num + 1)
            push!(chapters, (num=n, title=String(m.captures[2]), start_line=i))
        end
    end
end
println("Found $(length(chapters)) chapters:")
for c in chapters
    println("  $(c.num) $(c.title) @ line $(c.start_line)")
end

# ----------------------------------------------------------------------
# 2. Find chunk ranges (line spans of ```{julia} ... ``` blocks).
#    Each chunk: (open_line, close_line, content_lines, chapter_num)
# ----------------------------------------------------------------------
chunk_open_re  = r"^```\{julia\}"
chunk_close_re = r"^```$"

struct Chunk
    open_line::Int
    close_line::Int
    content::Vector{String}
    chapter::Int
end

function chapter_of(line::Int)
    chap = 0
    for c in chapters
        if c.start_line <= line
            chap = c.num
        end
    end
    return chap
end

chunks = Chunk[]
let i = 1
    while i <= n_lines
        if occursin(chunk_open_re, src_lines[i])
            j = i + 1
            while j <= n_lines && !occursin(chunk_close_re, src_lines[j])
                j += 1
            end
            content = src_lines[i+1:j-1]
            push!(chunks, Chunk(i, j, content, chapter_of(i)))
            i = j + 1
        else
            i += 1
        end
    end
end
println("\nFound $(length(chunks)) {julia} chunks")

# ----------------------------------------------------------------------
# 3. Preamble — the FIRST chunk (lines ~26-127): imports + compat shims.
#    Always copied verbatim to every chapter file.
# ----------------------------------------------------------------------
preamble_chunk = chunks[1]    # the imports/shims chunk
preamble_text  = join(preamble_chunk.content, "\n")
println("\nPreamble chunk: lines $(preamble_chunk.open_line)-$(preamble_chunk.close_line), $(length(preamble_chunk.content)) lines")

# ----------------------------------------------------------------------
# 4. Identify top-level helper definitions (function/const/struct/macro
#    definitions at indent 0) across all chunks. Skip the preamble
#    chunk — its contents (imports + shims) are always included.
#
#    For each helper:
#      - name (function/const/struct identifier)
#      - source text (with preceding docstring if any)
#      - chunk index where defined
#      - line range within the chunk's content (for excision)
# ----------------------------------------------------------------------
struct Helper
    name::String
    source::String           # docstring + definition
    chunk_idx::Int
    first_line_in_chunk::Int # 1-based within chunk.content
    last_line_in_chunk::Int
end

# Patterns for top-level definitions at indent 0:
const FUN_RE   = r"^function\s+([A-Za-z_][\w!]*)"     # function NAME(
# Single-line method: NAME(args) = expr. Disallow `"` inside args so we
# don't false-match e.g. `println("r(N)    = ", r_N)` where the string
# literal happens to contain `(...)= `. Also reject `==` and `=>`.
const ARROW_RE = r"^([A-Za-z_][\w!]*)\s*\([^)\"]*\)\s*=(?![=>])"
const CONST_RE = r"^const\s+([A-Za-z_][\w]*)"
const STRUCT_RE = r"^(?:mutable\s+)?struct\s+([A-Za-z_][\w]*)"
const MACRO_RE = r"^macro\s+([A-Za-z_][\w!]*)"

# Built-in names that should never be treated as user helpers even if a
# heuristic matches them (defensive — the regex tightening above should
# already exclude these, but belt-and-suspenders).
const BUILTIN_NAMES = Set([
    "println", "print", "printstyled", "@printf", "@sprintf",
    "show", "display", "error", "throw", "warn", "info",
    "return", "if", "for", "while", "let", "begin", "do",
    "true", "false", "nothing", "missing", "NaN", "Inf",
])

# Detect a docstring (triple-quoted string) starting at indent 0.
is_doc_open(line)  = startswith(line, "\"\"\"")
is_doc_close(line) = startswith(line, "\"\"\"")

helpers = Helper[]
helper_by_name = Dict{String, Helper}()

# Backfill: if the line just above a def is the closing of a """..."""
# docstring, capture the docstring as part of the def's source.
function backfill_docstring_start(lines, defidx)
    # Look one line above; skip exactly one blank if present.
    k = defidx - 1
    while k >= 1 && isempty(strip(lines[k]))
        k -= 1
    end
    k < 1 && return defidx
    # The line at k must be `"""` (closing) for a docstring above it.
    occursin(r"^\"\"\"\s*$", lines[k]) || return defidx
    # Walk back to find the matching opening `"""`.
    j = k - 1
    while j >= 1 && !occursin(r"^\"\"\"", lines[j])
        j -= 1
    end
    return j >= 1 ? j : defidx
end

for (ci, ch) in enumerate(chunks)
    ci == 1 && continue   # skip preamble
    lines = ch.content
    n = length(lines)
    idx = 1
    in_doc = false        # track whether we're inside a `"""..."""` block
    while idx <= n
        line = lines[idx]

        # Toggle docstring state on lines that contain `"""`. A line that
        # contains an *even* number of `"""` tokens doesn't change state
        # (e.g. a one-line `""" ... """` docstring). Skip every line while
        # inside a docstring so that prose like "function f(β) of ..."
        # is never matched as a real definition.
        n_triple = count("\"\"\"", line)
        if isodd(n_triple)
            in_doc = !in_doc
            idx += 1
            continue
        end
        if in_doc
            idx += 1
            continue
        end

        # Match a top-level definition at indent 0 (skip blanks/comments
        # implicitly — if none of the regexes match, we just advance).
        m_fun    = match(FUN_RE, line)
        m_arrow  = match(ARROW_RE, line)
        m_const  = match(CONST_RE, line)
        m_struct = match(STRUCT_RE, line)
        m_macro  = match(MACRO_RE, line)

        if m_fun !== nothing || m_struct !== nothing || m_macro !== nothing
            name = m_fun !== nothing ? m_fun.captures[1] :
                   m_struct !== nothing ? m_struct.captures[1] :
                   m_macro.captures[1]
            block_end = idx
            for k in idx+1:n
                if lines[k] == "end"
                    block_end = k
                    break
                end
            end
            if String(name) in BUILTIN_NAMES
                idx = block_end + 1
                continue
            end
            doc_start = backfill_docstring_start(lines, idx)
            src = join(lines[doc_start:block_end], "\n")
            h = Helper(String(name), src, ci, doc_start, block_end)
            push!(helpers, h)
            helper_by_name[h.name] = h
            idx = block_end + 1
        elseif m_arrow !== nothing
            name = m_arrow.captures[1]
            if String(name) in BUILTIN_NAMES
                idx += 1
                continue
            end
            block_end = idx
            if occursin(r"\b(begin|let|quote|do)\s*$", line)
                for k in idx+1:n
                    if lines[k] == "end"
                        block_end = k
                        break
                    end
                end
            end
            doc_start = backfill_docstring_start(lines, idx)
            src = join(lines[doc_start:block_end], "\n")
            h = Helper(String(name), src, ci, doc_start, block_end)
            push!(helpers, h)
            helper_by_name[h.name] = h
            idx = block_end + 1
        elseif m_const !== nothing
            name = m_const.captures[1]
            if String(name) in BUILTIN_NAMES
                idx += 1
                continue
            end
            doc_start = backfill_docstring_start(lines, idx)
            src = join(lines[doc_start:idx], "\n")
            h = Helper(String(name), src, ci, doc_start, idx)
            push!(helpers, h)
            helper_by_name[h.name] = h
            idx += 1
        else
            idx += 1
        end
    end
end
println("\nFound $(length(helpers)) top-level definitions (candidates, pre-filter)")
println("Sample: ", [h.name for h in first(helpers, 10)])

open("extracted_funcs_ch2.txt", "w") do io
    for h in helpers
        if h.chunk_idx == 2
            ch_open = chunks[2].open_line
            src_first = ch_open + h.first_line_in_chunk
            src_last  = ch_open + h.last_line_in_chunk
            println(io, h.name, "\t", src_first, "-", src_last, "\t", src_last - src_first + 1, " lines")
        end
    end
end
# Snapshot per-chunk candidate counts BEFORE filter
chunk_candidate_count = Dict{Int,Int}()
for h in helpers
    chunk_candidate_count[h.chunk_idx] = get(chunk_candidate_count, h.chunk_idx, 0) + 1
end
top_chunks = sort(collect(chunk_candidate_count), by = x->-x[2])[1:5]
println("Top chunks by helper candidates:")
for (ci, c) in top_chunks
    println("  chunk $ci  (lines $(chunks[ci].open_line)-$(chunks[ci].close_line))  $c helper candidates")
end

# ----------------------------------------------------------------------
# 5. Identify "content chunks" per chapter — chunks within a chapter's
#    line range that aren't dominated by helper definitions. A chunk is
#    a content chunk if it isn't ENTIRELY made up of top-level helpers
#    we already cataloged. (In practice a chunk is either a pure-helper
#    chunk or a pure-content chunk; mixed is rare.)
# ----------------------------------------------------------------------
function chunk_is_helper_only(ch::Chunk, ci::Int)
    # A chunk is helper-only if extracted helpers cover the vast majority
    # of meaningful lines (non-blank, non-comment). Two thresholds:
    #   - small chunks: uncovered <= 2 absolute (tolerate stragglers)
    #   - large chunks: uncovered / meaningful < 10% (giant Ch1 helper
    #     chunk is ~22k lines and may have a handful of stray @eval or
    #     `Base.parse(...)` lines that don't match the def heuristic but
    #     don't disqualify it from being a helper-only library chunk)
    hits = [h for h in helpers if h.chunk_idx == ci]
    isempty(hits) && return false
    covered_lines = Set{Int}()
    for h in hits
        for k in h.first_line_in_chunk:h.last_line_in_chunk
            push!(covered_lines, k)
        end
    end
    uncovered = 0
    meaningful = 0
    for (k, line) in enumerate(ch.content)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        meaningful += 1
        k in covered_lines && continue
        uncovered += 1
    end
    return uncovered <= 2 || uncovered / max(meaningful, 1) < 0.10
end

# Map chapter num -> Vector{Int} (indices of content chunks)
chapter_content_chunks = Dict{Int, Vector{Int}}()
for (ci, ch) in enumerate(chunks)
    ci == 1 && continue
    ch.chapter == 0 && continue
    chunk_is_helper_only(ch, ci) && continue
    push!(get!(chapter_content_chunks, ch.chapter, Int[]), ci)
end
for n in 1:18
    nc = length(get(chapter_content_chunks, n, Int[]))
    println("Chapter $(lpad(n,2)): $(nc) content chunks")
end

# ----------------------------------------------------------------------
# 5b. Two-pass filter: keep ONLY helpers that come from helper-only
#     chunks. Helpers defined inside content chunks (alongside calling
#     code) are chunk-local artifacts — they stay where they were and
#     should not be promoted to a global helpers block. This avoids
#     duplicating defs and also avoids label-clobbering across multiple
#     inline `commafmt`/`row` style chunk-local helpers.
# ----------------------------------------------------------------------
helper_only_chunk_set = Set{Int}()
for ci in eachindex(chunks)
    ci == 1 && continue
    if chunk_is_helper_only(chunks[ci], ci)
        push!(helper_only_chunk_set, ci)
    end
end
helpers = filter(h -> h.chunk_idx in helper_only_chunk_set, helpers)
empty!(helper_by_name)
for h in helpers
    helper_by_name[h.name] = h   # last def wins (rare; only matters for multi-method overloads)
end
println("After two-pass filter: $(length(helpers)) helpers from $(length(helper_only_chunk_set)) helper-only chunks")

# DIAGNOSTIC: per chunk, how many candidates were extracted and was it
# classified helper-only? (Pre-filter view, for debugging the heuristic.)
candidate_counts = Dict{Int, Int}()
for ci in eachindex(chunks)
    candidate_counts[ci] = 0
end
# Re-extract candidate count from the chunk-source-line span:
# count helpers in `pre_filter_helpers` would be cleaner, but we already
# overwrote `helpers`. Use a quick recount from the source positions.
# Show the top 10 candidate chunks by size:
chunk_def_counts = Dict{Int,Int}()
for ci in eachindex(chunks)
    ci == 1 && continue
    n = length(chunks[ci].content)
    chunk_def_counts[ci] = n
end
giant_chunks = sort(collect(chunk_def_counts), by = x->-x[2])[1:5]
println("\nLargest non-preamble chunks (by line count):")
for (ci, n) in giant_chunks
    ho = ci in helper_only_chunk_set
    n_helpers_here = count(h -> h.chunk_idx == ci, helpers)  # post-filter
    # Recompute counts via the function logic
    hits_pre = [h for h in helpers if h.chunk_idx == ci]
    cov = Set{Int}()
    for h in hits_pre
        for k in h.first_line_in_chunk:h.last_line_in_chunk
            push!(cov, k)
        end
    end
    uncov = 0; meaningful = 0
    uncov_lines_sample = String[]
    for (k, line) in enumerate(chunks[ci].content)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        meaningful += 1
        k in cov && continue
        uncov += 1
        length(uncov_lines_sample) < 10 && push!(uncov_lines_sample, "L$k: $(line)")
    end
    println("  chunk $ci  (lines $(chunks[ci].open_line)-$(chunks[ci].close_line))  $n total / $meaningful meaningful / $uncov uncovered / $(length(hits_pre)) defs  helper_only=$ho")
    if ci == 2
        println("  Sample uncovered lines in chunk 2:")
        for l in uncov_lines_sample
            println("    $l")
        end
    end
end

# ----------------------------------------------------------------------
# 6. Scan each helper's source for references to OTHER helpers
#    (dependency edges). Word-boundary match on helper names.
# ----------------------------------------------------------------------
function scan_refs(text::String, names::Vector{String})
    refs = Set{String}()
    for name in names
        # word-boundary match
        re = Regex("\\b" * name * "\\b")
        if occursin(re, text)
            push!(refs, name)
        end
    end
    return refs
end

all_names = unique!(sort([h.name for h in helpers]))
helper_deps = Dict{String, Set{String}}()
for h in helpers
    refs = scan_refs(h.source, all_names)
    delete!(refs, h.name)   # self-refs don't count
    helper_deps[h.name] = refs
end
println("\nDependency edges built for $(length(helper_deps)) helpers")

# Transitive closure helper
function closure(seeds::Set{String}, deps::Dict{String,Set{String}})
    out = Set{String}()
    queue = collect(seeds)
    while !isempty(queue)
        n = pop!(queue)
        n in out && continue
        push!(out, n)
        if haskey(deps, n)
            for d in deps[n]
                d in out || push!(queue, d)
            end
        end
    end
    return out
end

# ----------------------------------------------------------------------
# 7. For each chapter, compute the set of needed helpers, then write
#    chapter_NN.qmd with YAML + preamble + needed helpers + content.
# ----------------------------------------------------------------------
function chapter_safe_name(num::Int, title::String)
    cleaned = replace(title, r"[^\w\s-]" => "")
    cleaned = replace(cleaned, r"\s+" => "_")
    return "ch$(lpad(num,2,'0'))_$(cleaned).qmd"
end

# Get the YAML header from the source (lines 1..first ```{julia})
yaml_end = preamble_chunk.open_line - 1
yaml_lines = src_lines[1:yaml_end]
yaml_text  = join(yaml_lines, "\n")

# Build chapter -> seed helper names (referenced in content chunks)
for ch in chapters
    cn = ch.num
    haskey(chapter_content_chunks, cn) || continue
    content_idxs = chapter_content_chunks[cn]

    # Concatenate all content text for scanning
    content_text = join((join(chunks[i].content, "\n") for i in content_idxs), "\n\n")

    seeds = scan_refs(content_text, all_names)
    needed = closure(seeds, helper_deps)
    # Sort needed in their original definition order (preserve dependency order).
    # Exclude helpers whose source chunk is itself a preserved content chunk
    # of this chapter — those will appear naturally when we copy the chapter
    # content through (otherwise they'd be defined twice). Helpers in
    # helper-only chunks (mostly chunk 1's bulk helper section) are NOT in
    # `content_idxs`, so they always make it into the dedicated helpers chunk.
    preserved_for_ch = Set(content_idxs)
    ordered_helpers = [h for h in helpers if h.name in needed && !(h.chunk_idx in preserved_for_ch)]

    # Build chapter file content
    out_lines = String[]
    push!(out_lines, yaml_text)
    push!(out_lines, "")

    # Preamble {julia} chunk verbatim (imports + shims)
    push!(out_lines, "```{julia}")
    append!(out_lines, preamble_chunk.content)
    push!(out_lines, "```")
    push!(out_lines, "")

    # One {julia} chunk with the needed helpers (in definition order)
    if !isempty(ordered_helpers)
        push!(out_lines, "```{julia}")
        for (k, h) in enumerate(ordered_helpers)
            push!(out_lines, h.source)
            k < length(ordered_helpers) && push!(out_lines, "")
        end
        push!(out_lines, "```")
        push!(out_lines, "")
    end

    # Chapter header + content chunks (with original surrounding markdown
    # preserved by walking line ranges between chunk opens)
    # We walk the source from the chapter's start_line to the next
    # chapter's start_line and copy through markdown, but only KEEP
    # content chunks (skip helper chunks).
    cend = let nextchap = findfirst(c->c.num==cn+1, chapters)
        nextchap === nothing ? n_lines : chapters[nextchap].start_line - 1
    end
    cstart = ch.start_line
    # Build a set of helper-chunk index spans to skip
    helper_chunk_idxs = [ci for ci in eachindex(chunks)
                          if ci != 1 && chunks[ci].chapter == cn &&
                             chunk_is_helper_only(chunks[ci], ci)]
    skip_spans = [(chunks[ci].open_line, chunks[ci].close_line) for ci in helper_chunk_idxs]

    i = cstart
    while i <= cend
        # Are we inside a helper-chunk span to skip?
        in_skip = false
        for (a,b) in skip_spans
            if a <= i <= b
                i = b + 1
                in_skip = true
                break
            end
        end
        in_skip && continue
        push!(out_lines, src_lines[i])
        i += 1
    end

    fname = chapter_safe_name(cn, ch.title)
    open(fname, "w") do io
        write(io, join(out_lines, "\n"))
    end
    println("  wrote $(fname): $(length(ordered_helpers)) helpers, $(length(content_idxs)) content chunks")
end

println("\nDone.")
