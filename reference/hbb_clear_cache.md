# Clear Compiled Stan Model Cache

Removes all cached compiled Stan executables and hash files from the
persistent user-level cache directory. Also clears the in-memory session
cache so that the next call to a fitting function will trigger fresh
compilation.

## Usage

``` r
hbb_clear_cache()
```

## Value

Invisibly returns `TRUE` if the disk cache existed and was removed,
`FALSE` if no cache was found.

## See also

[`hbb_compile`](https://joonho112.github.io/hurdlebb/reference/hbb_compile.md)
to compile Stan models.

## Examples

``` r
if (FALSE) { # \dontrun{
hbb_clear_cache()
} # }
```
