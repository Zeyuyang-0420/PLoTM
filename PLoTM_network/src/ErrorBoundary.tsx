import { Component } from "react";
import type { ErrorInfo, ReactNode } from "react";

interface Props { children: ReactNode; onReset?: () => void }
interface State { error: Error | null }

/** Keeps a transient render error from blanking the whole network panel;
 *  shows a recover button instead of an empty iframe. */
export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };
  static getDerivedStateFromError(error: Error): State { return { error }; }
  componentDidCatch(error: Error, info: ErrorInfo) { console.error("network render error", error, info); }
  reset = () => { this.setState({ error: null }); this.props.onReset?.(); };
  render() {
    if (this.state.error) {
      return (
        <div className="graphwrap crash" role="alert">
          <div className="crashbox">
            <b>The graph hit a snag.</b>
            <p>{this.state.error.message}</p>
            <button onClick={this.reset}>recover</button>
          </div>
        </div>
      );
    }
    return this.props.children;
  }
}
