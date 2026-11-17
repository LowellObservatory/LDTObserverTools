# Repository Instructions

## Development Environment

- Use the `obstools` conda environment for development, tests, formatting, and
  other repository commands. For non-interactive commands, prefer
  `conda run -n obstools <command>` so the environment is explicit.

## Python Style

- Write full NumPy-style docstrings for every class, function, and method.
  Include all applicable sections, such as `Parameters`, `Returns`, `Yields`,
  `Raises`, `Attributes`, and `Notes`.
- Add type hints to every function and method, including return types.
- Minimize private helpers and inline short, single-use logic when that remains
  readable. Use a private helper when code is duplicated or when extracting a
  large block makes the calling routine's otherwise readable logic easier to
  follow.
- When creating a Python module, follow the repository's existing license-header
  convention. Use nearby modules as the source of truth for the MPL 2.0 header
  and preserve any applicable authorship or creation-date conventions.

## Efficient Use of Tokens

- Read only the files and sections needed for the current task; use `rg` to
  locate relevant definitions, tests, and conventions before opening files.
- Prefer concise progress updates and avoid repeating repository context,
  command output, or explanations already established in the conversation.
- Keep patches focused. Do not reformat or rewrite unrelated code, and inspect
  diffs rather than rereading entire files after small edits.
- Run the narrowest relevant tests first, expanding to broader checks only when
  the change or failures justify it.
- Summarize lengthy command output and report the actionable result rather than
  reproducing it verbatim.
- When a command fails, read the error before retrying. Change the retry only
  when the error suggests a specific correction; do not spend tokens repeating
  an identical call that is likely to fail for the same reason.
- Combine independent read-only inspections into one command when that keeps
  the output small and clear. Do not combine unrelated edits or validations
  whose failures would be difficult to diagnose.

## Sandbox Failures

- Treat sandbox errors as environment failures, not as evidence that the source
  code or test is broken. Preserve the original error for the final summary
  when it affects verification or completion.
- First check command syntax, the working directory, and whether the requested
  path is inside the writable workspace. Use explicit paths and avoid broad
  filesystem operations.
- If a necessary command fails because of sandboxing, permissions, networking,
  or package access, retry it once with the available approval or escalation
  mechanism and a narrowly scoped justification.
- Prefer a semantically equivalent, non-destructive tool when the preferred
  tool itself cannot start. Keep the fallback limited to the requested files
  and verify the resulting diff immediately.
- Do not work around the sandbox by weakening permissions, writing outside the
  approved workspace, exposing secrets, or using destructive commands.
- If both the approved retry and a safe fallback fail, stop retrying. Report the
  exact blocker, what was attempted, and what user action or environment change
  is required.
- After any fallback edit, inspect the changed file or diff before running tests
  so a no-op or partial edit is caught without wasting a test run.

## Git Commits

- When asked to commit changes, write the commit message as exactly three prose
  paragraphs: an imperative summary, a description of the implementation and
  user-visible behavior, and a summary of validation or other relevant context.
- After the three paragraphs, add a blank line followed by the trailer
  `Co-authored-by: Codex <codex@openai.com>`. Treat the trailer as commit
  metadata, not as one of the three prose paragraphs.
