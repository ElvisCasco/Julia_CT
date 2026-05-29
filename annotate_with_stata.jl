#!/usr/bin/env julia
# annotate_with_stata.jl
# For a given chapter qmd (default: ch01_Stata_basics.qmd), prepend the
# corresponding Stata source from `Cameron & Trivedi_stata.qmd` to each
# Julia chunk, wrapped as a `#= ... =#` reference block.
#
# Matching strategy: positional pairing within section. For each Julia
# chunk, walk up to the most recent markdown header — that section name
# is the bucket. Within each section, the n-th Julia chunk pairs with
# the n-th Stata chunk. If a section has fewer Stata chunks than Julia
# chunks (or vice versa), the surplus is left unannotated.

const STATA_QMD  = "Cameron & Trivedi_stata.qmd"
const TARGET_QMD = length(ARGS) >= 1 ? ARGS[1] : "ch01_Stata_basics.qmd"

# ----------------------------------------------------------------------
# Parse a qmd into chunks of language `lang`. Track the most recent
# markdown header as the "section" for each chunk.
# ----------------------------------------------------------------------
function parse_chunks(path::AbstractString, lang::AbstractString)
    src = readlines(path)
    chunks = NamedTuple{(:open, :close, :content, :section),
                        Tuple{Int,Int,Vector{String},String}}[]
    open_re  = Regex("^```\\{$(lang)\\}")
    close_re = r"^```\s*$"
    section  = ""
    let i = 1
        while i <= length(src)
            line = src[i]
            m = match(r"^(#+)\s+(.+)$", line)
            if m !== nothing
                section = strip(String(m.captures[2]))
                i += 1
            elseif occursin(open_re, line)
                j = i + 1
                while j <= length(src) && !occursin(close_re, src[j])
                    j += 1
                end
                push!(chunks, (open = i, close = j,
                               content = src[i+1:j-1], section = section))
                i = j + 1
            else
                i += 1
            end
        end
    end
    return chunks, src
end

stata_chunks, _   = parse_chunks(STATA_QMD,  "stata")
julia_chunks, jsrc = parse_chunks(TARGET_QMD, "julia")

println("Parsed $(length(stata_chunks)) Stata chunks from $(STATA_QMD)")
println("Parsed $(length(julia_chunks)) Julia chunks from $(TARGET_QMD)")

# ----------------------------------------------------------------------
# Group chunk indices by section, preserving order.
# ----------------------------------------------------------------------
function group_by_section(chs)
    d = Dict{String, Vector{Int}}()
    for (i, c) in enumerate(chs)
        push!(get!(d, c.section, Int[]), i)
    end
    return d
end

stata_by_section = group_by_section(stata_chunks)
julia_by_section = group_by_section(julia_chunks)

# ----------------------------------------------------------------------
# Pair Julia chunks with Stata chunks positionally within each section.
# ----------------------------------------------------------------------
pairs = Dict{Int, Int}()
unpaired_sections = String[]
for (sec, jidxs) in julia_by_section
    sidxs = get(stata_by_section, sec, Int[])
    if isempty(sidxs) && !isempty(jidxs)
        push!(unpaired_sections, sec)
    end
    for (k, jidx) in enumerate(jidxs)
        if k <= length(sidxs)
            pairs[jidx] = sidxs[k]
        end
    end
end

n_paired = length(pairs)
n_julia  = length(julia_chunks)
println("Paired $(n_paired)/$(n_julia) Julia chunks with Stata refs")
if !isempty(unpaired_sections)
    println("Sections with Julia chunks but no Stata chunks (sample):")
    for s in first(unpaired_sections, 8)
        println("  - $s")
    end
end

# ----------------------------------------------------------------------
# Rewrite TARGET_QMD, inserting `#= ... =#` Stata blocks right after each
# matched Julia chunk's opening ```{julia} fence.
# ----------------------------------------------------------------------
stata_for_open = Dict{Int, Vector{String}}()
for (jidx, sidx) in pairs
    stata_for_open[julia_chunks[jidx].open] = stata_chunks[sidx].content
end

new_lines = String[]
sizehint!(new_lines, length(jsrc) + 4 * n_paired)
for (i, line) in enumerate(jsrc)
    push!(new_lines, line)
    if haskey(stata_for_open, i)
        sc = stata_for_open[i]
        push!(new_lines, "#=")
        for s in sc
            push!(new_lines, s)
        end
        push!(new_lines, "=#")
    end
end

# Backup once
backup = TARGET_QMD * ".bak"
if !isfile(backup)
    open(backup, "w") do io
        write(io, join(jsrc, "\n"))
    end
    println("Backup saved → $backup")
end

open(TARGET_QMD, "w") do io
    write(io, join(new_lines, "\n"))
end
println("Updated $TARGET_QMD (+$(length(new_lines) - length(jsrc)) lines)")
