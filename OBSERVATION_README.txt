Artificial Universe — Observation (au watch)
============================================

6 files. Extract over your project and overwrite.
Also included, OUTSIDE au/:
  PROJECT_VISION_for_project_context.md
    -> upload to the project's Context panel, replacing PROJECT_VISION.md.
       Adds "Appearance Is Derived, Like Everything Else".

USAGE
  cargo build --release
  ./target/release/au watch --ticks 3000 --every 10 --out world.html
  ./target/release/au watch --config my.kv --ticks 5000 --out run.html
  open world.html   (self-contained; works offline, no JS, no CDN)

WHAT IT IS
  A read-only recorder + a one-file HTML/SVG view. Recorder::observe takes
  &World, so it CANNOT change the world - the borrow checker enforces it.
  Four panels: protocell census, populations (medium vs encapsulated, log),
  temperature band (signed log), matter across the boundary (gross).

THIS IS INSTRUMENTATION, NOT RENDERING
  Colours here are arbitrary and mean nothing. Real rendering must DERIVE
  appearance from what things are made of; nothing here does. Do not let the
  two blur.

FIRST FINDING
  3000 ticks -> 3000 nucleated, 531 lysed, 2 divided. The protocell world is
  a conveyor belt of newborns, not a reproducing population. Phase 5i step 2's
  headline test passes on 2 > 0; the claim holds, the magnitude is far weaker
  than its changelog implied. First item for Phase 5j.
